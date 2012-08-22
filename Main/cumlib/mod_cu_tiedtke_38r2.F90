module mod_cu_tiedtke_38r2

  use mod_constants
  use mod_dynparam

  private

  integer , parameter :: n_vmass = 0        ! Using or not vector mass
  real(dp) , parameter :: rlpal1 = 0.15D0   ! Smoothing coefficient
  real(dp) , parameter :: rlpal2 = 20.0D0   ! Smoothing coefficient
  real(dp) , parameter :: rmfsoluv = 1.0D0  ! Mass flux solver for momentum
  real(dp) , parameter :: rcucov = 0.05D0   ! Convective cloud cover for rain
                                            ! evporation
  real(dp) , parameter :: rcpecons = 5.44D-4/rgas
                                            ! Coefficient for rain evaporation
                                            ! below cloud
  real(dp) , parameter :: rtaumel = 5.0D0*3.6D3*1.5D0
                                            ! Relaxation time for melting of
                                            ! snow
  real(dp) , parameter :: rhebc = 0.8D0     ! Critical relative humidity below
                                            ! cloud at which evaporation starts
  real(dp) , parameter :: rmfcfl = 5.0D0    ! Massflux multiple of cfl stability
                                            ! criterium
  real(dp) , parameter :: detrpen = 0.75D-4 ! Detrainment rate for penetrative
                                            ! convection
  real(dp) , parameter :: entrorg = 1.75D-3 ! Entrainment for positively
                                            ! buoyant convection 1/(m)
  real(dp) , parameter :: entrdd = 3.0D-4   ! Average entrainment rate
                                            ! for downdrafts
  real(dp) , parameter :: rmfcmax = 1.0D0   ! Maximum massflux value allowed
                                            ! for updrafts etc
  real(dp) , parameter :: rmfcmin = 1.0D-10 ! Minimum massflux value (safety)
  real(dp) , parameter :: rmfdeps = 0.30D0  ! Fractional massflux for
                                            ! downdrafts at lfs
  real(dp) , parameter :: rdepths = 2.4D4   ! Maximum allowed cloud thickness
                                            ! for shallow cloud depth (Pa)
  real(dp) , parameter :: rvdifts = 1.5D0   ! Factor for time step weighting in
                                            ! *vdf....*
  real(dp) , parameter :: rmfsoltq = 1.0D0  ! Mass flux solver for T and q
  real(dp) , parameter :: zqmax = 0.5D0
  real(dp) , parameter :: rlmin = 1.0D-8

  logical , public :: lphylin = .false.  ! Linearized physics is activated ?
  logical , public :: lmfmid  ! True if midlevel convection is on
  logical , public :: lmfdd   ! True if cumulus downdraft is on
  logical , public :: lepcld  ! True if prognostic cloud scheme is on
  logical , public :: lmfdudv ! True if cumulus friction is on

  integer , parameter :: njkt1 = 2
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
!**** *CUSTRAT* - COMPUTES T,Q TENDENCIES FOR STRATOCUMULUS
!                 CONVECTION
!     M.TIEDTKE      E.C.M.W.F.    4/89 MODIF. 12/89
!     PURPOSE.
!     --------
!           THIS ROUTINE DOES THE PARAMETERIZATION OF BOUNDARY-LAYER
!     MIXING BY ENHANCED VERTICAL DIFFUSION OF SENSIBLE HEAT
!     AND MOISTURE FOR THE CASE OF STRATOCUMULUS CONVECTION.
!     THE SCHEME IS ONLY APPLIED IN THE BOUNDARY-LAYER AND ONLY
!     WHEN NEITHER PENETRATIVE NOR SHALLOW CONVECTION ARE ACTIVATED.
!**   INTERFACE.
!     ----------
!           THIS ROUTINE IS CALLED FROM *CUCALL*:
!     IT TAKES ITS INPUT FROM THE LONG-TERM STORAGE:
!     T,Q AT (T-1) AS WELL AS THE PROVISIONAL T,Q TENDENCIES AND
!     RETURNS ITS OUTPUT TO THE SAME SPACE:
!       MODIFIED TENDENCIES OF T AND Q
!     METHOD.
!     -------
!           ENHANCED VERTICAL DIFFUSION OF MOISTURE AND SENSIBLE
!     HEAT OCCURS, WHENEVER
!        1. LIFTED SURFACE AIR IS BUOYANT AND
!        2. CONDENSATION LEVEL EXISTS FOR FREE CONVECTION
!     THEN THE EXCHANGE COEFFICIENT IS AS FOLLOWS;
!         K=C1 IN CLOUD LAYER
!         K=C1*F(RH) AT CLOUD TOP (TOP ENTRAINMENT)
!     THE MATRIX INVERSION IS PERFORMED ANALOGOUSLY TO ROUTINE *VDIFF*
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    INPUT PARAMETERS (LOGICAL):
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
!C     INPUT PARAMETERS (REAL)
!    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS            PA
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PQSAT*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
!    *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S)
!    UPDATED PARAMETERS (REAL):
!    *PTENT*        TEMPERATURE TENDENCY                           K/S
!    *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)
!     EXTERNALS.
!     ----------
!          *CUADJTQ* ADJUST T AND Q DUE TO CONDENSATION IN ASCENT
!     Modifications.
!     --------------
!     G. Mozdzynski 2000-11-29: Corrections required for reproducibility
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!----------------------------------------------------------------------
!
  subroutine custrat(kidia,kfdia,klon,ktdia,klev,ldcum,ptsphy,pap,paph, &
                     pgeo,pten,pqen,pqsat,penth,ptent,ptenq)
    implicit none

    integer , intent(in) :: klon
    integer , intent(in) :: klev
    integer , intent(in) :: kidia
    integer , intent(in) :: kfdia
    integer :: ktdia ! argument not used
    logical , dimension(klon) , intent(in) :: ldcum
    real(dp) , intent(in) :: ptsphy
    real(dp) , dimension(klon,klev) , intent(in) :: pap , pgeo , pten , &
               pqen , pqsat
    real(dp) , dimension(klon,klev+1) , intent(in) :: paph
    real(dp) , dimension(klon,klev) , intent(out) :: penth
    real(dp) , dimension(klon,klev) , intent(inout) :: ptent , ptenq
    real(dp) , dimension(klon,klev) :: ztc , zqc , zcf , zcptgz , ztdif , &
               zqdif , zebs
    real(dp) , dimension(klon,klev+1) :: zap
    real(dp) , dimension(klon) :: zcpts , zqs , ztcoe , zqold , zpp
    integer , dimension(klon,klev) :: ilab
    logical , dimension(klon) :: llflag , llo2 , llbl
    integer :: icall , ik , ilevh , jk , jl
    real(dp) :: zbuo , zcons1 , zcons2 , zcons3 , zdisc , zdqdt , zdtdt , &
                zfac , zkdiff1 , zkdiff2 , zqdp , ztmst , ztpfac1 , ztpfac2
    !
    ! 1. Physical constants and parameters.
    !
    ztpfac1 = rvdifts
    ztpfac2 = d_one/ztpfac1
    zkdiff1 = 10.0D0
    zkdiff2 = 2.5D0
    ztmst = ptsphy
    zcons1 = ztpfac1*ztmst*egrav**d_two/(d_half*rgas)
    zcons2 = d_one/ztmst
    zcons3 = ztmst*cpd
    ilevh = klev/2
    !
    ! 2. Preliminary computations.
    !
    do jk = 1 , klev
      do jl = kidia , kfdia
        zcptgz(jl,jk) = pgeo(jl,jk)+pten(jl,jk)*cpd
        zcf(jl,jk) = d_zero
        ilab(jl,jk) = 0
      end do
    end do
    !
    ! 3. Determine exchange coefficients therefore
    !      (a) Lift surface air, check for buoyancy and set flag
    !      (b) Then define diffusion coefficients,i.e.
    !          k=c1 for cloud layer
    !          k=c1*f(rh) for cloud top (top entrainment)
    !
    do jl = kidia , kfdia
      ztc(jl,klev) = pten(jl,klev)+0.25D0
      zqc(jl,klev) = pqen(jl,klev)
      if ( .not. ldcum(jl) ) then
        ilab(jl,klev) = 1
      else
        ilab(jl,klev) = 0
      end if
      llo2(jl) = .false.
      llbl(jl) = .true.
    end do

    do jk = klev-1 , ilevh , -1
      do jl = kidia , kfdia
        if ( pap(jl,jk) < 0.9D0*paph(jl,klev+1) ) llbl(jl) = .false.
      end do
      do jl = kidia , kfdia
        if ( llbl(jl) ) then
          ztc(jl,jk) = (ztc(jl,jk+1)*cpd+pgeo(jl,jk+1)-pgeo(jl,jk))/cpd
          zqc(jl,jk) = zqc(jl,jk+1)
          if ( ilab(jl,jk+1) > 0 ) then
            llflag(jl) = .true.
          else
            llflag(jl) = .false.
          end if
          zap(jl,jk) = pap(jl,jk)
          zpp(jl) = pap(jl,jk)
          zqold(jl) = zqc(jl,jk)
        end if
      end do
      do jl = kidia , kfdia
        if ( .not.llbl(jl) ) llflag(jl) = .false.
      end do
      ik = jk
      icall = 1
      call cuadjtq(kidia,kfdia,klon,ktdia,klev,ik,zpp,ztc,zqc,llflag,icall)
      do jl = kidia , kfdia
        if ( llbl(jl) ) then
          if ( dabs(zqc(jl,jk)-zqold(jl)) > dlowval ) then
            ilab(jl,jk) = 2
          end if
        end if
      end do
!dir$ ivdep
!ocl novrec
      do jl = kidia , kfdia
        if ( llbl(jl) ) then
          zbuo = ztc(jl,jk)*(d_one+retv*zqc(jl,jk)) - &
                             pten(jl,jk)*(d_one+retv*pqen(jl,jk))
          if ( zbuo < d_zero ) ilab(jl,jk) = 0
          if ( zbuo > d_zero .and. &
               ilab(jl,jk) == 0 .and. ilab(jl,jk+1) == 1) ilab(jl,jk) = 1
          if ( ilab(jl,jk) == 2 ) llo2(jl) = .true.
        end if
      end do
    end do

    do jl = kidia , kfdia
      llbl(jl) = .true.
    end do

    do jk = klev-1 , ilevh , -1
      do jl = kidia , kfdia
        if ( pap(jl,jk) < 0.9D0*paph(jl,klev+1) ) llbl(jl) = .false.
      end do
      do jl = kidia , kfdia
        if ( llbl(jl) ) then
          if ( ilab(jl,jk) == 2 ) then
            zcf(jl,jk) = zkdiff1
            if ( ilab(jl,klev-2) == 0 ) then
              zcf(jl,jk) = zkdiff2
            end if
          else
            zcf(jl,jk) = d_zero
          end if
          if ( zcf(jl,jk+1) > d_zero .and. ilab(jl,jk) == 0 ) then
            zcf(jl,jk) = zcf(jl,jk+1)*5.0D0 *                             &
                         max(pqen(jl,jk+1)/pqsat(jl,jk+1)-0.8D0,d_zero) * &
                         max(pqen(jl,jk+1)/pqsat(jl,jk+1)-pqen(jl,jk) /   &
                         pqsat(jl,jk),d_zero)
            llbl(jl) = .false.
          end if
        end if
      end do
    end do
    !
    ! 4.7 Exchange coefficients.
    !
    do jk = ilevh , klev-1
      do jl = kidia , kfdia
        zcf(jl,jk) = zcf(jl,jk)*zcons1*paph(jl,jk+1) / &
                     ((pgeo(jl,jk)-pgeo(jl,jk+1)) *    &
                      (pten(jl,jk)+pten(jl,jk+1)))
      end do
    end do
    !
    ! 4.8 Dummy surface values of t and q at surface
    !
    do jl = kidia , kfdia
      zcpts(jl) = ztpfac2*zcptgz(jl,klev)
      zqs(jl) = ztpfac2*pqen(jl,klev)
    end do
    !
    ! 5. Solution of the vertical diffusion equation.
    !
    !  5.1 Setting of right hand sides.
    !
    do jk = ilevh , klev
      do jl = kidia , kfdia
        ztdif(jl,jk) = ztpfac2*zcptgz(jl,jk)
        zqdif(jl,jk) = ztpfac2*pqen(jl,jk)
      end do
    end do
    !
    !  5.2 Top layer elimination.
    !
    do jl = kidia , kfdia
      ztcoe(jl) = zcf(jl,ilevh)
      zqdp = d_one/(paph(jl,ilevh+1)-paph(jl,ilevh))
      zdisc = d_one/(d_one+zcf(jl,ilevh)*zqdp)
      zebs(jl,ilevh) = zdisc*(zcf(jl,ilevh)*zqdp)
      zqdif(jl,ilevh) = zdisc*zqdif(jl,ilevh)
      ztdif(jl,ilevh) = zdisc*ztdif(jl,ilevh)
    end do
    !
    !  5.3 Elimination for layers below
    !
    do jk = ilevh+1 , klev
      do jl = kidia , kfdia
        zqdp = d_one/(paph(jl,jk+1)-paph(jl,jk))
        zfac = ztcoe(jl)*zqdp
        ztcoe(jl) = zcf(jl,jk)
        zdisc = d_one/(d_one+zfac*(d_one-zebs(jl,jk-1))+zcf(jl,jk)*zqdp)
        zebs(jl,jk) = zdisc*(zcf(jl,jk)*zqdp)
        zqdif(jl,jk) = zdisc*(zqdif(jl,jk)+zfac*zqdif(jl,jk-1))
        ztdif(jl,jk) = zdisc*(ztdif(jl,jk)+zfac*ztdif(jl,jk-1))
      end do
    end do
    do jl = kidia , kfdia
      zqdif(jl,klev) = zqdif(jl,klev)+(zebs(jl,klev)*zqs(jl))
      ztdif(jl,klev) = ztdif(jl,klev)+(zebs(jl,klev)*zcpts(jl))
    end do
    !
    !  5.5 Back-substitution.
    !
    do jk = klev-1 , ilevh , -1
      do jl = kidia , kfdia
        zqdif(jl,jk) = zqdif(jl,jk)+(zebs(jl,jk)*zqdif(jl,jk+1))
        ztdif(jl,jk) = ztdif(jl,jk)+(zebs(jl,jk)*ztdif(jl,jk+1))
      end do
    end do
    !
    ! 6. Incrementation of t and q tendencies.
    !
    do jk = ilevh , klev
      do jl = kidia , kfdia
        zdqdt = (zqdif(jl,jk)-ztpfac2*pqen(jl,jk))*zcons2
        ptenq(jl,jk) = ptenq(jl,jk)+zdqdt
        zdtdt = (ztdif(jl,jk)-ztpfac2*zcptgz(jl,jk))/zcons3
        ptent(jl,jk) = ptent(jl,jk)+zdtdt
        penth(jl,jk) = (ztdif(jl,jk)-ztpfac2*zcptgz(jl,jk))*zcons2
      end do
    end do
  end subroutine custrat
!
!**** *CUDTDQ* - UPDATES T AND Q TENDENCIES, PRECIPITATION RATES
!                DOES GLOBAL DIAGNOSTICS
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!          P.BECHTOLD        E.C.M.W.F.     10/05
!**   INTERFACE.
!     ----------
!          *CUDTDQ* IS CALLED FROM *CUMASTR*
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
!    *KCTOP*        CLOUD TOP LEVEL
!    *KDTOP*        TOP LEVEL OF DOWNDRAFTS
!    INPUT PARAMETERS (LOGICAL):
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
!    *LDDRAF*       FLAG: .TRUE. FOR DOWNDRAFT LEVEL
!    INPUT PARAMETERS (REAL):
!    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PGEO*         GEOPOTENTIAL ON FULL LEVELS                   M2/S2
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
!    *PQSEN*        SATURATION ENV. SPEC. HUMIDITY (T+1)          KG/KG
!    *PLGLAC*       FLUX OF FROZEN CLOUDWATER IN UPDRAFTS         KG/(M2*S)
!    *PLUDE*        DETRAINED LIQUID WATER                        KG/(M3*S)
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS          J/(M2*S)
!    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS        J/(M2*S)
!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS            KG/(M2*S)
!    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS          KG/(M2*S)
!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS              KG/(M2*S)
!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP.                    KG/(M2*S)
!    *PDPMEL*       CHANGE IN PRECIP.-FLUXES DUE TO MELTING       KG/(M2*S)
!    UPDATED PARAMETERS (REAL):
!    *PTENT*        TEMPERATURE TENDENCY                           K/S
!    *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)
!    OUTPUT PARAMETERS (REAL):
!    *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S)
!----------------------------------------------------------------------
!               MODIFICATIONS
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!       96-09-20       : changed so can be used with diagnost
!                        cloud scheme      D.GREGORY
!       99-06-04       : Optimisation      D.SALMOND
!       03-08-28       : Clean up LINPHYS  P.BECHTOLD
!       05-10-13       : implicit solution P.BECHTOLD
!----------------------------------------------------------------------
!
  subroutine cudtdqn(kidia,kfdia,klon,ktdia,klev,ktopm2,ktype,kctop,kdtop, &
                     ldcum,lddraf,ptsphy,paph,pgeoh,pgeo,pten,ptenh,pqen,  &
                     pqenh,pqsen,plglac,plude,pmfu,pmfd,pmfus,pmfds,pmfuq, &
                     pmfdq,pmful,pdmfup,pdpmel,ptent,ptenq,penth)
    implicit none

    integer , intent(in) :: klon
    integer , intent(in) :: klev
    integer , intent(in) :: kidia
    integer , intent(in) :: kfdia
    integer :: ktdia ! argument not used
    integer , intent(in) :: ktopm2
    integer , dimension(klon) , intent(in) :: ktype
    integer , dimension(klon) , intent(in) :: kctop
    integer , dimension(klon) , intent(in) :: kdtop
    logical , dimension(klon) , intent(inout) :: ldcum
    logical , dimension(klon) , intent(in) :: lddraf
    real(dp) , intent(in)    :: ptsphy
    real(dp) , dimension(klon,klev+1) , intent(in) :: paph , pgeoh
    real(dp) , dimension(klon,klev) , intent(in) :: pgeo , pten , pqen , &
         ptenh , pqenh , pqsen , plglac , pmfu , pmfd , pmfus , pmfds ,  &
         pmfuq , pmfdq , pmful , pdmfup , pdpmel
    real(dp) , dimension(klon,klev) , intent(inout) :: plude , ptent , ptenq
    real(dp) , dimension(klon,klev) , intent(out) :: penth

    logical :: lltest
    integer :: jk , ik , jl
    real(dp) :: ztsphy , zorcpd , zalv , zoealfa , ztarg , &
                zzp , zgq , zgs , zgh , zs , zq
    real(dp) , dimension(klon,klev) :: zmfus , zmfuq , zmfds , zmfdq
    real(dp) , dimension(klon,klev) :: zdtdt , zdqdt , zdp , zb , zr1 , zr2
    logical , dimension(klon,klev) :: llcumbas

    real(dp) , parameter :: zimp = 1.0D0 - rmfsoltq

    !
    ! 1.0 Setup and initializations
    !
    ztsphy = d_one/ptsphy

    do jk = 1 , klev
      do jl = kidia , kfdia
        penth(jl,jk) = d_zero
      end do
    end do
    !
    ! mass-flux approach switched on for deep convection only
    ! in the tangent-linear and adjoint versions
    !
    do jl = kidia , kfdia
      if ( ktype(jl) /= 1 .and. lphylin ) ldcum(jl) = .false.
    end do
    !
    ! zero detrained liquid water if diagnostic cloud scheme to be used
    !
    ! this means that detrained liquid water will be evaporated in the
    ! cloud environment and not fed directly into a cloud liquid water
    ! variable

    !lltest = (.not.lepcld.and..not.lencld2).or.(lphylin.and..not.lencld2)
    lltest = .not. lepcld .or. lphylin

    if ( lltest ) then
      do jk = 1 , klev
        do jl = kidia , kfdia
          plude(jl,jk) = d_zero
        end do
      end do
    end if
    do jk = 1 , klev
      do jl = kidia , kfdia
        if ( ldcum(jl) ) then
          zdp(jl,jk) = egrav/(paph(jl,jk+1)-paph(jl,jk))
          zmfus(jl,jk) = pmfus(jl,jk)
          zmfds(jl,jk) = pmfds(jl,jk)
          zmfuq(jl,jk) = pmfuq(jl,jk)
          zmfdq(jl,jk) = pmfdq(jl,jk)
        end if
      end do
    end do
    !
    if ( rmfsoltq > d_zero ) then
      !
      ! 2.0 Recompute convective fluxes if implicit
      !
      do jk = ktopm2 , klev
        ik = jk-1
!dir$ ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
            !
            ! Compute interpolating coefficients zgs and zgq for
            ! half-level values
            !
            zgq = (pqenh(jl,jk)-pqen(jl,ik))/pqsen(jl,jk)
            zgh = cpd*pten(jl,jk)+pgeo(jl,jk)
            zgs = (cpd*(ptenh(jl,jk)-pten(jl,ik))+pgeoh(jl,jk)-pgeo(jl,ik))/zgh
            !
            ! half-level environmental values for s and q
            !
            zs = cpd*(zimp*pten(jl,ik)+zgs*pten(jl,jk)) + &
                 pgeo(jl,ik)+zgs*pgeo(jl,jk)
            zq = zimp*pqen(jl,ik)+zgq*pqsen(jl,jk)
            zmfus(jl,jk) = pmfus(jl,jk)-pmfu(jl,jk)*zs
            zmfuq(jl,jk) = pmfuq(jl,jk)-pmfu(jl,jk)*zq
            if ( lddraf(jl) .and. jk >= kdtop(jl) ) then
              zmfds(jl,jk) = pmfds(jl,jk)-pmfd(jl,jk)*zs
              zmfdq(jl,jk) = pmfdq(jl,jk)-pmfd(jl,jk)*zq
            end if
          end if
        end do
      end do
    end if
    !
    ! 3.0 Compute tendencies
    !
    do jk = ktopm2 , klev
      if ( jk < klev ) then
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            if ( lphylin ) then
              ztarg = pten(jl,jk)
              zoealfa = 0.545D0*(tanh(0.17D0*(ztarg-mpcrt))+d_one)
              zalv = zoealfa*wlhv+(d_one-zoealfa)*wlhs
            else
              zalv = foelhmcu(pten(jl,jk))
            end if
            zdtdt(jl,jk) = zdp(jl,jk)*rcpd*(zmfus(jl,jk+1)-zmfus(jl,jk) + &
                                            zmfds(jl,jk+1)-zmfds(jl,jk) + &
                                  wlhf*plglac(jl,jk)-wlhf*pdpmel(jl,jk) - &
                                      zalv*(pmful(jl,jk+1)-pmful(jl,jk) - &
                                            plude(jl,jk)-pdmfup(jl,jk)))
            zdqdt(jl,jk) = zdp(jl,jk)*(zmfuq(jl,jk+1)-zmfuq(jl,jk) + &
                                       zmfdq(jl,jk+1)-zmfdq(jl,jk) + &
                                       pmful(jl,jk+1)-pmful(jl,jk) - &
                                       plude(jl,jk)-pdmfup(jl,jk))
          end if
        end do
      else
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            if ( lphylin ) then
              ztarg = pten(jl,jk)
              zoealfa = 0.545D0*(tanh(0.17D0*(ztarg-mpcrt))+d_zero)
              zalv = zoealfa*wlhv+(d_one-zoealfa)*wlhs
            else
              zalv = foelhmcu(pten(jl,jk))
            end if
            zdtdt(jl,jk) = -zdp(jl,jk)*rcpd*(zmfus(jl,jk)+zmfds(jl,jk) + &
                         wlhf*pdpmel(jl,jk)-zalv*(pmful(jl,jk)+pdmfup(jl,jk)))
            zdqdt(jl,jk) = -zdp(jl,jk)*(zmfuq(jl,jk)+zmfdq(jl,jk) + &
                                       (pmful(jl,jk)+pdmfup(jl,jk)))
          end if
        end do
      end if
    end do
    if ( dabs(rmfsoltq) < dlowval ) then
      !
      ! 3.1 Update tendencies
      !
      do jk = ktopm2 , klev
        do jl = kidia , kfdia
          if( ldcum(jl) ) then
            ptent(jl,jk) = ptent(jl,jk)+zdtdt(jl,jk)
            ptenq(jl,jk) = ptenq(jl,jk)+zdqdt(jl,jk)
            penth(jl,jk) = zdtdt(jl,jk)*cpd
          end if
        end do
      end do
    else
      !
      ! 3.2 Implicit solution
      !
      ! fill bi-diagonal matrix vectors a=k-1, b=k, c=k+1;
      ! reuse zmfus=a
      ! zdtdt and zdqdt correspond to the rhs ("constants") of the equation
      ! the solution is in zr1 and zr2
      llcumbas(:,:) = .false.
      zb(:,:) = d_one
      zmfus(:,:) = d_zero
      !
      ! fill vectors a, b and rhs
      !
      do jk = ktopm2 , klev
        ik = jk+1
        do jl = kidia , kfdia
          llcumbas(jl,jk) = ldcum(jl) .and. jk >= kctop(jl)-1
          if ( llcumbas(jl,jk) ) then
            zzp = rmfsoltq*zdp(jl,jk)*ptsphy
            zmfus(jl,jk) = -zzp*(pmfu(jl,jk)+pmfd(jl,jk))
            zdtdt(jl,jk) = zdtdt(jl,jk)*ptsphy+pten(jl,jk)
            zdqdt(jl,jk) = zdqdt(jl,jk)*ptsphy+pqen(jl,jk)
            ! zdtdt(jl,jk) = (zdtdt(jl,jk)+ptent(jl,jk))*ptsphy+pten(jl,jk)
            ! zdqdt(jl,jk) = (zdqdt(jl,jk)+ptenq(jl,jk))*ptsphy+pqen(jl,jk)
            if ( jk < klev ) then
              zb(jl,jk) = d_one+zzp*(pmfu(jl,ik)+pmfd(jl,ik))
            else
              zb(jl,jk) = d_one
            end if
          end if
        end do
      end do

      call cubidiag(kidia,kfdia,klon,klev,kctop,llcumbas,zmfus,zb,zdtdt,zr1)
      call cubidiag(kidia,kfdia,klon,klev,kctop,llcumbas,zmfus,zb,zdqdt,zr2)
      !
      ! Compute tendencies
      !
      do jk = ktopm2 , klev
        do jl = kidia , kfdia
          if ( llcumbas(jl,jk) ) then
            ptent(jl,jk) = ptent(jl,jk)+(zr1(jl,jk)-pten(jl,jk))*ztsphy
            ptenq(jl,jk) = ptenq(jl,jk)+(zr2(jl,jk)-pqen(jl,jk))*ztsphy
            ! ptent(jl,jk)=(zr1(jl,jk)-pten(jl,jk))*ztsphy
            ! ptenq(jl,jk)=(zr2(jl,jk)-pqen(jl,jk))*ztsphy
            penth(jl,jk) = (zr1(jl,jk)-pten(jl,jk))*ztsphy
          end if
        end do
      end do
    end if
  end subroutine cudtdqn
!
!          M.TIEDTKE         E.C.M.W.F.     7/86 MODIF. 12/89
!          PURPOSE
!          -------
!          THIS ROUTINE DOES THE FINAL CALCULATION OF CONVECTIVE
!          FLUXES IN THE CLOUD LAYER AND IN THE SUBCLOUD LAYER
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
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
!    *KDTOP*        TOP LEVEL OF DOWNDRAFTS
!    INPUT PARAMETERS (LOGICAL):
!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
!    INPUT PARAMETERS (REAL):
!    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS            PA
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    UPDATED PARAMETERS (INTEGER):
!    *KTYPE*        SET TO ZERO IF LDCUM=.FALSE.
!    UPDATED PARAMETERS (LOGICAL):
!    *LDDRAF*       SET TO .FALSE. IF LDCUM=.FALSE. OR KDTOP<KCTOP
!    UPDATED PARAMETERS (REAL):
!    *PMFU*         MASSFLUX IN UPDRAFTS                          KG/(M2*S)
!    *PMFD*         MASSFLUX IN DOWNDRAFTS                        KG/(M2*S)
!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS          J/(M2*S)
!    *PMFDS*        FLUX OF DRY STATIC ENERGY IN DOWNDRAFTS        J/(M2*S)
!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS            KG/(M2*S)
!    *PMFDQ*        FLUX OF SPEC. HUMIDITY IN DOWNDRAFTS          KG/(M2*S)
!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS              KG/(M2*S)
!    *PLUDE*        DETRAINED LIQUID WATER                        KG/(M3*S)
!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS        KG/(M2*S)
!    *PDMFDP*       FLUX DIFFERENCE OF PRECIP. IN DOWNDRAFTS      KG/(M2*S)
!    *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                    KG/(M2*S)
!    OUTPUT PARAMETERS (REAL):
!    *PDPMEL*       CHANGE IN PRECIP.-FLUXES DUE TO MELTING       KG/(M2*S)
!    *PLGLAC*       FLUX OF FROZEN CLOUD WATER IN UPDRAFTS        KG/(M2*S)
!    *PMFLXR*       CONVECTIVE RAIN FLUX                          KG/(M2*S)
!    *PMFLXS*       CONVECTIVE SNOW FLUX                          KG/(M2*S)
!    *PRAIN*        TOTAL PRECIP. PRODUCED IN CONV. UPDRAFTS      KG/(M2*S)
!                   (NO EVAPORATION IN DOWNDRAFTS)
!          EXTERNALS
!          ---------
!          NONE
!          MODIFICATIONS
!          -------------
!             99-06-14 : Optimisation        D.SALMOND
!             03-08-28 : Clean up LINPHYS    P.BECHTOLD
!                        Bugs in Evapor.
!             05-02-11 : Extend DDMflux to   P.BECHTOLD
!                        surface if zero
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!----------------------------------------------------------------------
!
  subroutine cuflxn(kidia,kfdia,klon,ktdia,klev,ptsphy,pten,pqen,pqsen,  &
                    ptenh,pqenh,paph,pap,pgeoh,ldland,ldcum,kcbot,kctop, &
                    kdtop,ktopm2,ktype,lddraf,pmfu,pmfd,pmfus,pmfds,     &
                    pmfuq,pmfdq,pmful,plude,pdmfup,pdmfdp,pdpmel,plglac, &
                    pmflxr,pmflxs,prain,pmfdde_rate)
    implicit none

    integer , intent(in) :: klon
    integer , intent(in) :: klev
    integer , intent(in) :: kidia
    integer , intent(in) :: kfdia
    integer :: ktdia ! argument not used
    real(dp) , intent(in)    :: ptsphy
    real(dp) , dimension(klon,klev) , intent(in) :: pten , pqen , ptenh , &
             pqenh , pap
    real(dp) , dimension(klon,klev+1) , intent(in) :: paph , pgeoh
    real(dp) , dimension(klon,klev) , intent(inout) :: pqsen
    logical , dimension(klon) :: ldland ! argument not used
    logical , dimension(klon) , intent(in) :: ldcum
    integer , dimension(klon) , intent(in) :: kcbot , kctop , kdtop
    integer , intent(out) :: ktopm2
    integer , dimension(klon) , intent(inout) :: ktype
    logical , dimension(klon) , intent(inout) :: lddraf
    real(dp) , dimension(klon,klev) , intent(inout) :: pmfu , pmfd , &
           pmfus , pmfds , pmfuq , pmfdq , pmful , pdmfup , pdmfdp , &
           plglac , pmfdde_rate
    real(dp) , dimension(klon,klev) , intent(out) :: plude , pdpmel
    real(dp) , dimension(klon,klev+1) , intent(out) :: pmflxr , pmflxs
    real(dp) , dimension(klon) , intent(out) :: prain

    real(dp) , dimension(klon) :: zrhebc
    integer :: ik , ikb , jk , jl
    integer , dimension(klon) :: idbas
    logical :: llddraf

    real(dp) :: zalfaw , zcons1 , zcons1a , zcons2 , zdenom , zdrfl ,  &
       zdrfl1 , zfac , zfoeewi , zfoeewl , zoealfa , zoeewm , zoelhm , &
       zpdr , zpds , zrfl , zrfln , zrmin , zrnew , zsnmlt , ztarg ,   &
       ztmst , zzp
    !
    ! Specify constants
    !
    ztmst = ptsphy
    zcons1a = cpd/(wlhf*egrav*rtaumel)
    ! zcons2 = d_one/(egrav*ztmst)
    zcons2 = rmfcfl/(egrav*ztmst)
    !
    ! 1.0 Determine final convective fluxes
    !
    do jl = kidia , kfdia
      prain(jl) = d_zero
      if ( .not. ldcum(jl) .or. kdtop(jl) < kctop(jl) ) lddraf(jl) = .false.
      if ( .not. ldcum(jl) ) ktype(jl) = 0
      idbas(jl) = klev
      if ( ldland(jl) ) then
        ! zrhebc(jl) = rhebc
        zrhebc(jl) = 0.7D0
      else
        ! zrhebc(jl) = rhebc
        zrhebc(jl) = 0.9D0
      end if
    end do
    !
    ! To get identical results for different nproma force ktopm2 to 2
    !
    ktopm2 = 2
    do jk = ktopm2 , klev
!dir$ ivdep
!ocl novrec
      ikb = min(jk+1,klev)
      do jl = kidia , kfdia
        pmflxr(jl,jk) = d_zero
        pmflxs(jl,jk) = d_zero
        pdpmel(jl,jk) = d_zero
        if ( ldcum(jl) .and. jk >= kctop(jl) ) then
          pmfus(jl,jk) = pmfus(jl,jk)-pmfu(jl,jk) * &
                         (cpd*ptenh(jl,jk)+pgeoh(jl,jk))
          pmfuq(jl,jk) = pmfuq(jl,jk)-pmfu(jl,jk)*pqenh(jl,jk)
          plglac(jl,jk) = pmfu(jl,jk)*plglac(jl,jk)
          llddraf = lddraf(jl) .and. jk >= kdtop(jl)
          if ( llddraf ) then
            pmfds(jl,jk) = pmfds(jl,jk)-pmfd(jl,jk) * &
                           (cpd*ptenh(jl,jk)+pgeoh(jl,jk))
            pmfdq(jl,jk) = pmfdq(jl,jk)-pmfd(jl,jk)*pqenh(jl,jk)
          else
            pmfd(jl,jk) = d_zero
            pmfds(jl,jk) = d_zero
            pmfdq(jl,jk) = d_zero
            pdmfdp(jl,jk-1) = d_zero
          end if
          if ( llddraf .and. pmfd(jl,jk) < d_zero .and. &
               dabs(pmfd(jl,ikb)) < dlowval ) then
            idbas(jl) = jk
          end if
        else
          pmfu(jl,jk) = d_zero
          pmfd(jl,jk) = d_zero
          pmfus(jl,jk) = d_zero
          pmfds(jl,jk) = d_zero
          pmfuq(jl,jk) = d_zero
          pmfdq(jl,jk) = d_zero
          pmful(jl,jk) = d_zero
          plglac(jl,jk) = d_zero
          pdmfup(jl,jk-1) = d_zero
          pdmfdp(jl,jk-1) = d_zero
          plude(jl,jk-1) = d_zero
        end if
      end do
    end do
    pmflxr(:,klev+1) = d_zero
    pmflxs(:,klev+1) = d_zero
    !
    ! 1.5 Scale fluxes below cloud base linear dcrease
    !
!dir$ ivdep
!ocl novrec
    do jl = kidia , kfdia
      if ( ldcum(jl) ) then
        ikb = kcbot(jl)
        ik = ikb+1
        zzp = ((paph(jl,klev+1)-paph(jl,ik))/(paph(jl,klev+1)-paph(jl,ikb)))
        if ( ktype(jl) == 3 ) zzp = zzp*zzp
        pmfu(jl,ik) = pmfu(jl,ikb)*zzp
        if ( lphylin ) then
          ztarg = ptenh(jl,ikb)
          zoealfa = 0.545D0*(tanh(0.17D0*(ztarg-mpcrt))+d_one)
          zoelhm = zoealfa*wlhv+(d_one-zoealfa)*wlhs
          pmfus(jl,ik) = (pmfus(jl,ikb)-zoelhm*pmful(jl,ikb))*zzp
        else
          pmfus(jl,ik) = (pmfus(jl,ikb)-foelhmcu(ptenh(jl,ikb)) * &
                          pmful(jl,ikb))*zzp
        end if
        pmfuq(jl,ik) = (pmfuq(jl,ikb)+pmful(jl,ikb))*zzp
        pmful(jl,ik) = d_zero
      end if
    end do
    do jk = ktopm2 , klev
!dir$ ivdep
!ocl novrec
      do jl = kidia , kfdia
        if ( ldcum(jl) .and. jk > kcbot(jl)+1 ) then
          ikb = kcbot(jl)+1
          zzp = ((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
          if ( ktype(jl) == 3 ) zzp = zzp*zzp
          pmfu(jl,jk) = pmfu(jl,ikb)*zzp
          pmfus(jl,jk) = pmfus(jl,ikb)*zzp
          pmfuq(jl,jk) = pmfuq(jl,ikb)*zzp
          pmful(jl,jk) = d_zero
        end if
        ik = idbas(jl)
        llddraf = lddraf(jl) .and. jk > ik .and. ik < klev
        if ( llddraf .and. ik == kcbot(jl)+1 ) then
          zzp = ((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ik)))
          if ( ktype(jl) == 3 ) zzp = zzp*zzp
          pmfd(jl,jk) = pmfd(jl,ik)*zzp
          pmfds(jl,jk) = pmfds(jl,ik)*zzp
          pmfdq(jl,jk) = pmfdq(jl,ik)*zzp
          pmfdde_rate(jl,jk) = -(pmfd(jl,jk-1)-pmfd(jl,jk))
        else if ( llddraf .and. ik /= kcbot(jl)+1 .and. jk == ik+1 ) then
          pmfdde_rate(jl,jk) = -(pmfd(jl,jk-1)-pmfd(jl,jk))
        end if
      end do
    end do
    !
    ! 2. Calculate rain/snow fall rates
    !    Calculate melting of snow
    !    Calculate evaporation of precip
    !
    do jk = ktopm2 , klev
      do jl = kidia , kfdia
        if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
          prain(jl) = prain(jl)+pdmfup(jl,jk)
          if ( pmflxs(jl,jk) > d_zero .and. pten(jl,jk) > tzero ) then
            zcons1 = zcons1a*(d_one+d_half*(pten(jl,jk)-tzero))
            zfac = zcons1*(paph(jl,jk+1)-paph(jl,jk))
            zsnmlt = min(pmflxs(jl,jk),zfac*(pten(jl,jk)-tzero))
            pdpmel(jl,jk) = zsnmlt
            if ( lphylin ) then
              ztarg=pten(jl,jk)-zsnmlt/zfac
              zoealfa = 0.545D0*(tanh(0.17D0*(ztarg-mpcrt))+d_one)
              zfoeewl = c2es*exp(c3les*(ztarg-tzero)/(ztarg-c4les))
              zfoeewi = c2es*exp(c3ies*(ztarg-tzero)/(ztarg-c4ies))
              zoeewm = zoealfa*zfoeewl + (d_one-zoealfa)*zfoeewi
              pqsen(jl,jk) = zoeewm/pap(jl,jk)
            else
              pqsen(jl,jk) = foeewmcu(pten(jl,jk)-zsnmlt/zfac)/pap(jl,jk)
            end if
          end if
          if ( lphylin ) then
            zalfaw = 0.545D0*(tanh(0.17D0*(pten(jl,jk)-mpcrt))+d_one)
          else
            zalfaw=foealfcu(pten(jl,jk))
          end if
          !
          ! No liquid precipitation above melting level
          !
          if ( pten(jl,jk) < tzero .and. zalfaw > d_zero ) then
            plglac(jl,jk) = plglac(jl,jk)+zalfaw*(pdmfup(jl,jk)+pdmfdp(jl,jk))
            zalfaw = d_zero
          end if
          pmflxr(jl,jk+1) = pmflxr(jl,jk)+zalfaw * &
                            (pdmfup(jl,jk)+pdmfdp(jl,jk))+pdpmel(jl,jk)
          pmflxs(jl,jk+1) = pmflxs(jl,jk)+(d_one-zalfaw) * &
                            (pdmfup(jl,jk)+pdmfdp(jl,jk))-pdpmel(jl,jk)
          if ( pmflxr(jl,jk+1)+pmflxs(jl,jk+1) < d_zero ) then
            pdmfdp(jl,jk) = -(pmflxr(jl,jk)+pmflxs(jl,jk)+pdmfup(jl,jk))
            pmflxr(jl,jk+1) = d_zero
            pmflxs(jl,jk+1) = d_zero
            pdpmel(jl,jk)   = d_zero
          else if ( pmflxr(jl,jk+1) < d_zero ) then
            pmflxs(jl,jk+1) = pmflxs(jl,jk+1)+pmflxr(jl,jk+1)
            pmflxr(jl,jk+1) = d_zero
          else if ( pmflxs(jl,jk+1) < d_zero ) then
            pmflxr(jl,jk+1) = pmflxr(jl,jk+1)+pmflxs(jl,jk+1)
            pmflxs(jl,jk+1) = d_zero
          end if
        end if
      end do
    end do
    !
    ! Reminder for conservation:
    ! pdmfup(jl,jk)+pdmfdp(jl,jk) = pmflxr(jl,jk+1)+pmflxs(jl,jk+1) - &
    !                               pmflxr(jl,jk)  -pmflxs(jl,jk)
    do jk = ktopm2 , klev
      do jl = kidia , kfdia
        if ( ldcum(jl) .and. jk >= kcbot(jl) ) then
          zrfl = pmflxr(jl,jk)+pmflxs(jl,jk)
          if ( zrfl > dlowval ) then
            zdrfl1 = rcpecons*max(d_zero,pqsen(jl,jk)-pqen(jl,jk))*rcucov * &
                  (sqrt(paph(jl,jk)/paph(jl,klev+1))/5.09D-3 * &
                   zrfl/rcucov)**0.5777D0*(paph(jl,jk+1)-paph(jl,jk))
            zrnew = zrfl-zdrfl1
            zrmin = zrfl-rcucov*max(d_zero,zrhebc(jl) * &
                pqsen(jl,jk)-pqen(jl,jk))*zcons2*(paph(jl,jk+1)-paph(jl,jk))
            zrnew = max(zrnew,zrmin)
            zrfln = max(zrnew,d_zero)
            zdrfl = min(d_zero,zrfln-zrfl)
            if ( lphylin ) then
              zalfaw = 0.545D0*(tanh(0.17D0*(pten(jl,jk)-mpcrt))+d_one)
            else
              zalfaw=foealfcu(pten(jl,jk))
            end if
            if ( pten(jl,jk) < tzero ) zalfaw = d_zero
            zpdr = zalfaw*pdmfdp(jl,jk)
            zpds = (d_one-zalfaw)*pdmfdp(jl,jk)
            zdenom = d_one/max(dlowval,pmflxr(jl,jk)+pmflxs(jl,jk))
            pmflxr(jl,jk+1) = pmflxr(jl,jk)+zpdr + &
                              pdpmel(jl,jk)+zdrfl*pmflxr(jl,jk)*zdenom
            pmflxs(jl,jk+1) = pmflxs(jl,jk)+zpds - &
                              pdpmel(jl,jk)+zdrfl*pmflxs(jl,jk)*zdenom
            pdmfup(jl,jk) = pdmfup(jl,jk)+zdrfl
            if ( pmflxr(jl,jk+1)+pmflxs(jl,jk+1) < d_zero ) then
              pdmfup(jl,jk) = pdmfup(jl,jk)-(pmflxr(jl,jk+1)+pmflxs(jl,jk+1))
              pmflxr(jl,jk+1) = d_zero
              pmflxs(jl,jk+1) = d_zero
              pdpmel(jl,jk)   = d_zero
            else if ( pmflxr(jl,jk+1) < d_zero ) then
              pmflxs(jl,jk+1) = pmflxs(jl,jk+1)+pmflxr(jl,jk+1)
              pmflxr(jl,jk+1) = d_zero
            else if ( pmflxs(jl,jk+1) < d_zero ) then
              pmflxr(jl,jk+1) = pmflxr(jl,jk+1)+pmflxs(jl,jk+1)
              pmflxs(jl,jk+1) = d_zero
            end if
          else
            pmflxr(jl,jk+1) = d_zero
            pmflxs(jl,jk+1) = d_zero
            pdmfdp(jl,jk)   = d_zero
            pdpmel(jl,jk)   = d_zero
          end if
        end if
      end do
    end do
  end subroutine cuflxn
!
!          THIS ROUTINE CALCULATES CLOUD BASE FIELDS
!          CLOUD BASE HEIGHT AND CLOUD TOP HEIGHT
!          A. Pier Siebesma   KNMI ********      
!          modified C Jakob (ECMWF) (01/2001) 
!          modified P Bechtold (ECMWF) (08/2002) 
!          (include cycling over levels to find unstable departure/base level+
!           mixed layer properties +w Trigger)
!          PURPOSE.
!          --------
!          TO PRODUCE CLOUD BASE AND CLOUD TOP VALUES FOR CU-PARAMETRIZATION
!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUMASTR*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD FIELDS VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=0 FOR STABLE LAYERS
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CLOUD LEVELS LEVEL
!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD TOP
!          (ENTRAINING PLUME, WITH ENTRAINMENT PROPORTIONAL TO (1/Z))
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    INPUT PARAMETERS (REAL):
! not used at the moment because we want to use linear intepolation
! for fields on the half levels.
!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS           K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PQSEN*        PROVISIONAL ENVIRONMENT SATU. HUMIDITY (T+1)  KG/KG
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    UPDATED PARAMETERS (REAL):
!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PUU*          U-VELOCITY IN UPDRAFTS                         M/S
!    *PVU*          V-VELOCITY IN UPDRAFTS                         M/S
!    UPDATED PARAMETERS (INTEGER):
!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CLOUD LEVELS
!    OUTPUT PARAMETERS (LOGICAL):
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDSC*         FLAG: .TRUE. IF BL-CLOUDS EXIST
!    OUTPUT PARAMETERS (INTEGER):
!    *KCBOT*       CLOUD BASE LEVEL !    
!    *KCTOP*       CLOUD TOP LEVEL = HEIGHEST HALF LEVEL 
!                  WITH A NON-ZERO CLOUD UPDRAFT.
!    *KBOTSC*      CLOUD BASE LEVEL OF BL-CLOUDS
!    *KDPL*        DEPARTURE LEVEL
!    *PCAPE*       PSEUDOADIABATIQUE max CAPE (J/KG)
!          EXTERNALS
!          ---------
!          *CUADJTQ* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             02-11-02 : Use fixed last possible departure level and 
!                        last updraft computation level for bit-reproducibility
!                                            D.Salmond &  J. Hague
!             03-07-03 : Tuning for p690     J. Hague
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!----------------------------------------------------------------------
!
  subroutine cubasen(kidia,kfdia,klon,ktdia,klev,ldland,ptenh,pqenh,   &
                     pgeoh,paph,pqhfl,pahfs,pten,pqen,pqsen,pgeo,puen, &
                     pven,ptu,pqu,plu,puu,pvu,pwubase,klab,ldcum,ldsc, &
                     kcbot,kbotsc,kctop,kdpl,pcape)
    implicit none

    integer , intent(in) :: klon 
    integer , intent(in) :: klev 
    integer , intent(in) :: kidia 
    integer , intent(in) :: kfdia 
    integer :: ktdia ! argument not used
    logical , dimension(klon) , intent(in) :: ldland
    real(dp) , dimension(klon,klev) , intent(in) :: ptenh , pqenh , &
              pten , pqen , pqsen , pgeo , puen , pven
    real(dp) , dimension(klon,klev+1) , intent(in) :: pgeoh , paph , &
              pqhfl , pahfs
    real(dp) , dimension(klon,klev) , intent(inout) :: ptu , pqu , &
              plu , puu , pvu
    real(dp) , dimension(klon) , intent(out) :: pwubase , pcape
    integer , dimension(klon,klev) , intent(inout) :: klab
    logical , dimension(klon) , intent(inout) :: ldcum
    logical , dimension(klon) , intent(out) :: ldsc
    integer , dimension(klon) , intent(inout) :: kcbot
    integer , dimension(klon) , intent(out) :: kbotsc , kctop , kdpl

    integer , dimension(klon) ::  ictop , icbot , ibotsc , idpl
    integer , dimension(klon,klev) :: ilab
    logical , dimension(klon) :: ll_ldbase , llgo_on , lldeep , lldcum , &
             lldsc , llfirst , llresetjl
    logical :: llreset
    integer :: icall , ik , ikb , is , jk , jl , jkk , jkt1 , jkt2 , jkt , jkb
    real(dp) , dimension(klon,klev) :: zs , zsuh , zwu2h , zbuoh , zlu , &
             zqu , ztu , zuu , zvu , zcape
    real(dp) , dimension(klon,klev+1) :: zsenh , zqenh
    real(dp) ,dimension(klon) :: zqold , zph , zmix , zdz , zcbase , &
             ztven1 , ztvu1 , zdtvtrig
    real(dp) :: zrho      ! density at surface (kg/m^3) 
    real(dp) :: zkhvfl    ! surface buoyancy flux (k m/s)
    real(dp) :: zws       ! sigma_w at lowest model halflevel (m/s)
    real(dp) :: zqexc     ! humidity excess at lowest model halflevel (kg/kg)
    real(dp) :: ztexc     ! temperature excess at lowest model halflevel (k)
    real(dp) :: zeps      ! fractional entrainment rate   [m^-1]
    real(dp) :: ztvenh    ! environment virtual temperature at half levels (k)  
    real(dp) :: ztvuh     ! updraft virtual temperature at half levels     (k)
    real(dp) :: zlglac    ! updraft liquid water frozen in one layer
    real(dp) :: zqsu , zcor , zdq , zalfaw , zfacw , zfaci , zfac ,  &
                zesdp , zdqsdt , zdtdp , zdp ,zpdifftop ,zpdiffbot , &
                zsf , zqf , zbuof , zz , ztmp
    real(dp) :: ztven2 , ztvu2 ! pseudoadiabatique t_v
    real(dp) :: zwork1 , zwork2 ! work arrays for t and w perturbations

    real(dp) , parameter :: zc2 = 0.55D0
    real(dp) , parameter :: zaw = 1.0D0
    real(dp) , parameter :: zbw = 1.0D0
    real(dp) , parameter :: zepsadd = 1.D-4
    !
    ! 0. Initialize fields
    !
    do jl = kidia , kfdia
      pwubase(jl) = d_zero
      llgo_on(jl) = .true.
      llfirst(jl) = .true.
      kdpl(jl) = klev
    end do

    jkt1 = njkt1
    jkt2 = njkt2

    do jk = 1 , klev
      do jl = kidia , kfdia
        ztu(jl,jk) = ptu(jl,jk)
        zqu(jl,jk) = pqu(jl,jk)
        zlu(jl,jk) = plu(jl,jk)
        zuu(jl,jk) = puu(jl,jk)
        zvu(jl,jk) = pvu(jl,jk)
        ilab(jl,jk) = klab(jl,jk)
        zcape(jl,jk) = d_zero
      end do
    end do
    !
    ! 1.1 Prepare fields on half levels by linear interpolation
    !     of specific humidity and static energy
    !
    do jk = 1 , klev
      do jl = kidia , kfdia
        zwu2h(jl,jk) = d_zero
        zs(jl,jk) = cpd*pten(jl,jk) + pgeo(jl,jk)
        zqenh(jl,jk) = pqenh(jl,jk)
        zsenh(jl,jk) = cpd*ptenh(jl,jk)+pgeoh(jl,jk)
      end do
    end do
    do jkk = klev , jkt1 , -1 ! big external loop for level testing:
      !
      ! find first departure level that produces deepest cloud top
      ! or take surface level for shallow convection and sc
      ! 1.2  Initialise fields at departure half model level
      !
      is = 0
      do jl = kidia , kfdia
        if ( llgo_on(jl) ) then
          is = is+1
          idpl(jl) = jkk      ! departure level
          icbot(jl) = jkk     ! cloud base level for conv., (-1 if not found)
          ibotsc(jl) = klev-1 ! base level for sc-clouds  , (-1 if not found)
          ictop(jl)  = klev-1 ! cloud top for convection (-1 if not found)
          lldcum(jl) = .false.    ! on exit: true if cloudbase=found
          lldsc (jl) = .false.    ! on exit: true if cloudbase=found
          ll_ldbase(jl) = .false. ! on exit: true if cloudbase=found
          zdtvtrig(jl) = d_zero
          zuu(jl,jkk) = puen(jl,jkk)*(paph(jl,jkk+1)-paph(jl,jkk))
          zvu(jl,jkk) = pven(jl,jkk)*(paph(jl,jkk+1)-paph(jl,jkk))
        end if 
      end do

      if ( is /= 0 ) then
        if ( jkk == klev ) then
          do jl = kidia , kfdia
            if ( llgo_on(jl) ) then
              zrho = paph(jl,jkk+1)/(rgas*(pten(jl,jkk) * &
                     (d_one+retv*pqen(jl,jkk))))
              zkhvfl = (pahfs(jl,jkk+1)*rcpd + &
                       retv*pten(jl,jkk)*pqhfl(jl,jkk+1))/zrho
              zws = 0.001D0 - 1.5D0*vonkar*zkhvfl * &
                    (pgeoh(jl,klev)-pgeoh(jl,klev+1))/pten(jl,klev)
              if ( zkhvfl < d_zero ) then
                zws = 1.2D0*zws**.3333D0
                ilab(jl,jkk) = 1
                ztexc = max(-1.5D0*pahfs(jl,jkk+1)/(zrho*zws*cpd),d_zero)
                zqexc = max(-1.5D0*pqhfl(jl,jkk+1)/(zrho*zws),d_zero)
                zqu(jl,jkk) = zqenh(jl,jkk) + zqexc
                zsuh(jl,jkk) = zsenh(jl,jkk) + cpd*ztexc
                ztu(jl,jkk) = (zsenh(jl,jkk)-pgeoh(jl,jkk))*rcpd + ztexc
                zlu(jl,jkk) = d_zero
                zwu2h(jl,jkk) = zws**d_two
                !
                !  Determine buoyancy at lowest half level
                !
                ztvenh = (d_one+retv*zqenh(jl,jkk)) * &
                         (zsenh(jl,jkk)-pgeoh(jl,jkk))*rcpd  
                ztvuh = (d_one+retv*zqu(jl,jkk))*ztu(jl,jkk)
                zbuoh(jl,jkk) = (ztvuh-ztvenh)*egrav/ztvenh
              else
                llgo_on(jl) = .false.      ! non-convective point
              end if
            end if
          end do
        else
          do jl = kidia , kfdia
            if ( llgo_on(jl) ) then
              zrho = paph(jl,jkk+1)/(rgas*(pten(jl,jkk) * &
                     (d_one+retv*pqen(jl,jkk))))
              ilab(jl,jkk) = 1
              ztexc = 0.2D0
              zqexc = 1.D-4
              zqu(jl,jkk) = zqenh(jl,jkk) + zqexc
              zsuh(jl,jkk) = zsenh(jl,jkk) + cpd*ztexc
              ztu(jl,jkk) = (zsenh(jl,jkk)-pgeoh(jl,jkk))*rcpd + ztexc
              zlu(jl,jkk) = d_zero
              ! construct mixed layer for parcels emanating in lowest 60 hpa
              if ( paph(jl,klev+1)-paph(jl,jkk-1) < 60.D2 ) then
                zqu(jl,jkk)  = d_zero
                zsuh(jl,jkk) = d_zero
                zwork1       = d_zero
                do jk = jkk+1 , jkk-1 , -1
                  if ( zwork1 < 50.D2 ) then
                    zwork2 = paph(jl,jk)-paph(jl,jk-1)
                    zwork1 = zwork1+zwork2
                    zqu(jl,jkk) = zqu(jl,jkk) + zqenh(jl,jk)*zwork2
                    zsuh(jl,jkk) = zsuh(jl,jkk)+zsenh(jl,jk)*zwork2
                  end if
                end do
                zqu(jl,jkk) = zqu(jl,jkk) / zwork1+zqexc
                zsuh(jl,jkk) = zsuh(jl,jkk)/zwork1+cpd*ztexc
                ztu(jl,jkk) = (zsuh(jl,jkk)-pgeoh(jl,jkk))*rcpd+ztexc
              end if
              zwu2h(jl,jkk) = d_one
              !
              !  determine buoyancy at lowest half level
              !
              ztvenh = (d_one+retv*zqenh(jl,jkk)) * &
                       (zsenh(jl,jkk)-pgeoh(jl,jkk))*rcpd  
              ztvuh = (d_one+retv*zqu(jl,jkk))*ztu(jl,jkk)
              zbuoh(jl,jkk) = (ztvuh-ztvenh)*egrav/ztvenh
            end if
          end do
        end if
      end if
      !
      ! 2.0 Do ascent in subcloud and layer,
      !     check for existence of condensation level,
      !     adjust t,q and l accordingly in *cuadjtq*,
      !     check for buoyancy and set flags
      !
      !
      ! 1.2  Do the vertical ascent until velocity becomes negative
      !
      do jk = jkk-1 , jkt2 , -1
        is = 0
        if ( jkk == klev ) then ! 1/z mixing for shallow
          do jl = kidia , kfdia
            if ( llgo_on(jl) ) then
              is = is+1
              zdz(jl) = (pgeoh(jl,jk) - pgeoh(jl,jk+1))*regrav
              zeps = zc2/((pgeoh(jl,jk)-pgeoh(jl,klev+1))*regrav) + zepsadd
              zmix(jl) = d_half*zdz(jl)*zeps
              zqf = (pqenh(jl,jk+1) + pqenh(jl,jk))*d_half
              zsf = (zsenh(jl,jk+1) + zsenh(jl,jk))*d_half
              ztmp = d_one/(d_one+zmix(jl))
              zqu(jl,jk) = (zqu(jl,jk+1)*(d_one-zmix(jl)) + &
                           d_two*zmix(jl)*zqf) * ztmp  
              zsuh(jl,jk) = (zsuh(jl,jk+1)*(d_one-zmix(jl)) + &
                            d_two*zmix(jl)*zsf) * ztmp  
              zqold(jl) = zqu(jl,jk)
              ztu(jl,jk) = (zsuh(jl,jk)-pgeoh(jl,jk))*rcpd
              zph(jl) = paph(jl,jk)
            end if
          end do
        else
          do jl = kidia , kfdia
            if ( llgo_on(jl) ) then
              is = is+1
              zdz(jl) = (pgeoh(jl,jk) - pgeoh(jl,jk+1))*regrav
              zqf = (pqenh(jl,jk+1) + pqenh(jl,jk))*d_half
              zsf = (zsenh(jl,jk+1) + zsenh(jl,jk))*d_half
              zmix(jl) = 0.4D0*entrorg*zdz(jl) * &
                         min(d_one,(pqsen(jl,jk)/pqsen(jl,klev))**d_three)
              zqu(jl,jk) = zqu(jl,jk+1)*(d_one-zmix(jl))+ zqf*zmix(jl)
              zsuh(jl,jk) = zsuh(jl,jk+1)*(d_one-zmix(jl))+ zsf*zmix(jl)
              zqold(jl) = zqu(jl,jk)
              ztu(jl,jk) = (zsuh(jl,jk)-pgeoh(jl,jk))*rcpd
              zph(jl) = paph(jl,jk)
            end if
          end do
        end if
        if ( is == 0 ) exit
        ik = jk
        icall = 1
        call cuadjtq(kidia,kfdia,klon,ktdia,klev,ik,zph,ztu,zqu,llgo_on,icall)
!dir$ ivdep
!ocl novrec
        do jl = kidia , kfdia
          if ( llgo_on(jl) ) then
            !
            ! Add condensation to water
            !
            zdq = max(zqold(jl)-zqu(jl,jk),d_zero)
            zlu(jl,jk) = zlu(jl,jk+1)+zdq
            !
            ! Freezing
            !
            zlglac = zdq*((d_one-foealfcu(ztu(jl,jk))) - &
                          (d_one-foealfcu(ztu(jl,jk+1))))  
            !
            ! Pseudo-microphysics
            !
            if ( jkk==klev ) then  ! no precip for shallow
              zlu(jl,jk) = min(zlu(jl,jk),5.D-3)
              !
              ! Chose a more pseudo-adiabatic formulation as original
              ! overestimates water loading effect and therefore strongly
              ! underestimates cloud thickness
              !
            else 
              zlu(jl,jk) = d_half*zlu(jl,jk) 
            end if
            !
            ! Update dry static energy after condensation + freezing
            !
            zsuh(jl,jk) = cpd*(ztu(jl,jk)+wlhfocp*zlglac)+pgeoh(jl,jk)
            !
            ! Buoyancy on half and full levels
            !
            ztvuh = (d_one+retv*zqu(jl,jk)-zlu(jl,jk)) * &
                    ztu(jl,jk)+wlhfocp*zlglac  
            ztvenh = (d_one+retv*zqenh(jl,jk)) * &
                     (zsenh(jl,jk)-pgeoh(jl,jk))*rcpd  
            zbuoh(jl,jk) = (ztvuh-ztvenh)*egrav/ztvenh
            zbuof = (zbuoh(jl,jk) + zbuoh(jl,jk+1))*d_half
            !
            ! Solve kinetic energy equation
            !
            ztmp = d_one/(d_one+d_two*zbw*zmix(jl))
            zwu2h(jl,jk) = (zwu2h(jl,jk+1)*(d_one-d_two*zbw*zmix(jl)) + &
                           d_two*zaw*zbuof*zdz(jl)) * ztmp  
            !
            ! Compute pseudoadiabatique cape for diagnostics
            !
            ztvu2 = ztu(jl,jk)   *(d_one+retv*zqu(jl,jk))
            ztven2 = ptenh(jl,jk)*(d_one+retv*pqenh(jl,jk))
            if ( jk == jkk-1 ) then
              ztvu1(jl) = ztvu2
              ztven1(jl) = ztven2
            end if
            zbuof = (ztvu2+ztvu1(jl)-ztven1(jl)-ztven2)/ztven2
            zbuof = zbuof*zdz(jl)*egrav
            zcape(jl,jkk) = zcape(jl,jkk) + max(d_zero,zbuof)
            ztvu1(jl) = ztvu2
            ztven1(jl) = ztven2
            !
            ! First layer with liquid water - find exact cloud base
            !
            if ( zlu(jl,jk) > d_zero .and. ilab(jl,jk+1) == 1 ) then
              ik = jk+1
              zqsu = foeewm(ztu(jl,ik))/paph(jl,ik)
              zqsu = min(d_half,zqsu)
              zcor = d_one/(d_one-retv*zqsu)
              zqsu = zqsu*zcor
              zdq = min(d_zero,zqu(jl,ik)-zqsu)
              zalfaw = foealfa(ztu(jl,ik))
              zfacw = c5les/((ztu(jl,ik)-c4les)**d_two)
              zfaci = c5ies/((ztu(jl,ik)-c4ies)**d_two)
              zfac = zalfaw*zfacw+(d_one-zalfaw)*zfaci
              zesdp = foeewm(ztu(jl,ik))/paph(jl,ik)
              zcor = d_one/(d_one-retv*zesdp)
              zdqsdt = zfac*zcor*zqsu
              zdtdp = rgas*ztu(jl,ik)/(cpd*paph(jl,ik))
              zdp = zdq/(zdqsdt*zdtdp)
              zcbase(jl) = paph(jl,ik)+zdp
              !
              ! Chose nearest half level as cloud base
              !
              zpdifftop = zcbase(jl)-paph(jl,jk)
              zpdiffbot = paph(jl,jk+1)-zcbase(jl)
              if ( zpdifftop > zpdiffbot .and. &
                   zwu2h(jl,jk+1) > d_zero ) then
                jkb = min(klev-1,jk+1)
                ilab(jl,jkb) = 2 
                ilab(jl,jk)  = 2
                ll_ldbase(jl) = .true.
                lldsc(jl)     = .true.
                ibotsc(jl) = jkb
                icbot(jl)  = jkb
                zlu(jl,jk+1) = rlmin
              else if ( zpdifftop <= zpdiffbot .and. &
                        zwu2h(jl,jk) > d_zero) then
                ilab(jl,jk) = 2
                ll_ldbase(jl) = .true.
                lldsc(jl)     = .true.
                ibotsc(jl) = jk
                icbot(jl)  = jk
              end if
              jkb = icbot(jl)
            end if
            !
            ! Decide on presence of convection, cloud base and
            ! cloud top based on kinetic energy
            !
            if ( zwu2h(jl,jk) < d_zero ) then
              llgo_on(jl) = .false.             
              if ( zlu(jl,jk+1) > d_zero ) then
                ictop(jl) = jk
                lldcum(jl) = .true.
              else
                lldcum(jl) = .false.
              end if
            else
              if ( zlu(jl,jk) > d_zero ) then
                ilab(jl,jk) = 2
              else
                ilab(jl,jk) = 1
              end if
            end if
          end if
        end do
        if ( lmfdudv .and. jkk == klev ) then
          do jl = kidia , kfdia
            if ( .not. ll_ldbase(jl) .and. llgo_on(jl) ) then
              zuu(jl,jkk) = zuu(jl,jkk)+puen(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
              zvu(jl,jkk) = zvu(jl,jkk)+pven(jl,jk)*(paph(jl,jk+1)-paph(jl,jk))
            end if
          end do
        end if
        ! if ( is == 0 ) exit
      end do
      if ( jkk == klev ) then
        !
        ! Set values for departure level for pbl clouds = first model level
        !
        do jl = kidia , kfdia
          ldsc(jl) = lldsc(jl)
          if ( ldsc(jl) ) then
            kbotsc(jl) = ibotsc(jl)
          else
            kbotsc(jl) = -1
          end if
          llgo_on(jl) = .false.
          jkt = ictop(jl)
          jkb = icbot(jl)
          lldeep(jl) = paph(jl,jkb)-paph(jl,jkt) > rdepths
          if ( lldeep(jl) ) lldcum(jl) = .false. ! no deep allowed for klev
          lldeep(jl) = .false. ! for deep convection start only at level klev-1
                               ! and form mixed layer, so go on
          ! test further for deep convective columns as not yet found
          if ( lldeep(jl) ) llfirst(jl) = .false.
          llgo_on(jl) = .not. lldeep(jl)
          if ( lldcum(jl) ) then
            kcbot(jl) = icbot(jl)
            kctop(jl) = ictop(jl)
            kdpl(jl)  = idpl(jl)
            ldcum(jl) = lldcum(jl)
            pwubase(jl) = sqrt(max(zwu2h(jl,jkb),d_zero))
          else
            kctop(jl) = -1
            kcbot(jl) = -1
            kdpl(jl) = klev-1
            ldcum(jl) = .false.
            pwubase(jl) = d_zero
          end if
        end do
        do jk = klev , 1 , -1
          do jl = kidia , kfdia
            jkt = ictop(jl)
            if ( jk >= jkt ) then
              klab(jl,jk) = ilab(jl,jk)
              ptu(jl,jk) = ztu(jl,jk)
              pqu(jl,jk) = zqu(jl,jk)
              plu(jl,jk) = zlu(jl,jk)
            end if
          end do
        end do
      end if
      if ( jkk < klev ) then
        llreset = .false.
        do jl = kidia , kfdia
          if ( .not. lldeep(jl) ) then
            jkt = ictop(jl)
            jkb = icbot(jl)
            ! test on cloud thickness and buoyancy
            lldeep(jl) = paph(jl,jkb)-paph(jl,jkt) >= rdepths 
            ! lldeep(jl)=paph(jl,jkb)-paph(jl,jkt)>=rdepths &
            !            .and. zdtvtrig(jl)>0._jprb
          end if
          llresetjl(jl) = lldeep(jl) .and. llfirst(jl)
          llreset = llreset .or. llresetjl(jl)
        end do
        if ( llreset ) then
          do jk = klev , 1 , -1
            do jl = kidia , kfdia
              ! keep first departure level that produces deep cloud
              ! if ( lldeep(jl) .and. llfirst(jl) ) then 
              if ( llresetjl(jl) ) then 
                jkt = ictop(jl)
                jkb = idpl(jl)
                if ( jk <= jkb .and. jk >= jkt ) then
                  klab(jl,jk) = ilab(jl,jk)
                  ptu(jl,jk) = ztu(jl,jk)
                  pqu(jl,jk) = zqu(jl,jk)
                  plu(jl,jk) = zlu(jl,jk)
                else 
                  klab(jl,jk) = 1
                  ptu(jl,jk) = ptenh(jl,jk)
                  pqu(jl,jk) = pqenh(jl,jk)
                  plu(jl,jk) = d_zero
                end if
                if ( jk < jkt ) klab(jl,jk) = 0
              end if
            end do
          end do
        end if
        do jl = kidia , kfdia
          if ( lldeep(jl) .and. llfirst(jl) ) then
            kdpl(jl)   = idpl(jl)
            kctop(jl)  = ictop(jl)
            kcbot(jl)  = icbot(jl)
            ldcum(jl)  = lldcum(jl)
            ldsc(jl)   = .false.
            kbotsc(jl) = -1
            jkb = kcbot(jl)
            pwubase(jl) = sqrt(max(zwu2h(jl,jkb),d_zero))
            ! No initialization of wind for deep here, this is done in
            ! cuini and cuascn
            llfirst(jl) = .false.
          end if
          llgo_on(jl) = .not. lldeep(jl)
        end do
      end if
    end do ! end of big loop for search of departure level     
    !
    ! Chose maximum cape value
    !
    do jl = kidia , kfdia
      pcape(jl) = maxval(zcape(jl,:))
    end do
  end subroutine cubasen
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
