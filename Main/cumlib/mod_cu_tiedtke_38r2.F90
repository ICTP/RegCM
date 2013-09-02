module mod_cu_tiedtke_38r2

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam

  private

  public :: sucumf , custrat , cumastrn , cuancape2

  integer(ik4) , parameter :: n_vmass = 0   ! Using or not vector mass
  real(rk8) , parameter :: rlpal1 = 0.15D0  ! Smoothing coefficient
  real(rk8) , parameter :: rlpal2 = 20.0D0  ! Smoothing coefficient
  real(rk8) , parameter :: rcucov = 0.05D0  ! Convective cloud cover for rain
                                            ! evporation
  real(rk8) , parameter :: rcpecons = 5.44D-4/rgas
                                          ! Coefficient for rain evaporation
                                          ! below cloud
  real(rk8) :: rtau
  real(rk8) , parameter :: rtaumel = 5.0D0*3.6D3*1.5D0
                                          ! Relaxation time for melting of
                                          ! snow
  real(rk8) , parameter :: rhebc = 0.8D0  ! Critical relative humidity below
                                          ! cloud at which evaporation starts
  real(rk8) :: rmfcfl ! Massflux multiple of cfl stability criterium
  real(rk8) , parameter :: rmflic = 1.0D0 ! Use CFL mass flux limit (1) or
                                          ! absolut limit (0)
  real(rk8) , parameter :: rmflia = 0.0D0 ! Value of absolut mass flux limit
  real(rk8) , parameter :: rmfsoluv = 1.0D0 ! Mass flux solver for momentum
  real(rk8) , parameter :: rmfsoltq = 1.0D0 ! Mass flux solver for T and q
  real(rk8) , parameter :: rmfsolct = 1.0D0 ! Mass flux solver for chem tracers
  real(rk8) , parameter :: ruvper = 0.3D0 ! Updraught velocity perturbation for
                                          ! implicit (m/s)
  real(rk8) , parameter :: rprcon = 1.4D03 ! coefficients for determining
                                           ! conversion from cloud water
  real(rk8) , parameter :: detrpen = 0.75D-4 ! Detrainment rate for penetrative
                                             ! convection
  real(rk8) , parameter :: entrorg = 1.75D-3 ! Entrainment for positively
                                             ! buoyant convection 1/(m)
  real(rk8) , parameter :: entshalp = 2.0D0  ! shallow entrainment factor for
                                             ! entrorg
  real(rk8) , parameter :: entrdd = 3.0D-4   ! Average entrainment rate
                                             ! for downdrafts
  real(rk8) , parameter :: rmfcmax = 1.0D0   ! Maximum massflux value allowed
                                             ! for updrafts etc
  real(rk8) , parameter :: rmfcmin = 1.0D-10 ! Minimum massflux value (safety)
  real(rk8) , parameter :: rmfdeps = 0.30D0  ! Fractional massflux for
                                             ! downdrafts at lfs
  real(rk8) , parameter :: rdepths = 2.4D4   ! Maximum allowed cloud thickness
                                             ! for shallow cloud depth (Pa)
  real(rk8) , parameter :: rvdifts = 1.5D0   ! Factor for time step weighting in
                                             ! *vdf....*
  real(rk8) , parameter :: zqmax = 0.5D0
  real(rk8) , parameter :: rlmin = 1.0D-8

  logical , public :: lmfmid = .true.   ! True if midlevel convection is on
  logical , public :: lmfdd = .true.    ! True if cumulus downdraft is on
  logical , public :: lepcld = .true.   ! True if prognostic cloud scheme is on
  logical , public :: lmfdudv = .true.  ! True if cumulus friction is on
  logical , public :: lmfscv = .true.   ! True if shallow convection is on
  logical , public :: lmfpen = .true.   ! True if penetrative convection is on
  logical , public :: lmfuvdis = .true. ! use kinetic energy dissipation
                                        ! (addit T-tendency)
  logical , public :: lmftrac = .true.  ! Convective chemical tracer transport
  logical , public :: lmfsmooth = .false. ! Smoothing of mass fluxes 
                                          ! top/bottom for tracers
  logical , public :: lmfwstar = .false. ! Grant w* closure for shallow conv.

  integer(ik4) :: njkt1 , njkt2 , njkt3 , njkt4 , njkt5 , njkt6

  contains
!
!     THIS ROUTINE DEFINES DISPOSABLE PARAMETERS FOR MASSFLUX SCHEME
!
!          M.TIEDTKE         E.C.M.W.F.    2/89
!
!          INTERFACE
!          ---------
!
!          THIS ROUTINE IS CALLED FROM *INIPHY*
!
!          MODIFICATIONS
!          -------------
!
  subroutine sucumf(ksmax,klev,pmean)
    implicit none
    integer(ik4) , intent(in) :: ksmax , klev
    real(rk8) , dimension(klev) , intent(in) :: pmean
    integer(ik4) :: nflevg , jlev
    nflevg = klev
    rtau = d_one+264.0D0/real(ksmax)
    rtau = min(3.0D0,rtau)
    if ( ksmax >= 511 ) then
      rmfcfl = 3.0D0
    else
      rmfcfl = 5.0D0
    end if
    njkt1 = 2
    njkt2 = 2
    njkt3 = nflevg-2
    do jlev = nflevg , 2 , -1
      if ( pmean(jlev)/pmean(klev)*stdp > 350.D2 ) njkt1 = jlev
      if ( pmean(jlev)/pmean(klev)*stdp >  60.D2 ) njkt2 = jlev
      if ( pmean(jlev)/pmean(klev)*stdp > 950.D2 ) njkt3 = jlev
      if ( pmean(jlev)/pmean(klev)*stdp > 850.D2 ) njkt4 = jlev
      if ( pmean(jlev)/pmean(klev)*stdp > 500.D2 ) njkt5 = jlev
      if ( pmean(jlev)/pmean(klev)*stdp > 700.D2 ) njkt6 = jlev
    end do
    njkt3 = min(nflevg-2,njkt3)
  end subroutine sucumf
!
!***** CUANCAPE2 - COMPUTE APPROXIMATE CAPE,CIN  USING THETAE AND THETAES
!
!     E. HOLM + P. BECHTOLD     E.C.M.W.F.     13/10/2005
!
!     PURPOSE 
!     -------
!                 ESTIMATE CAPE FIRST FOR A MIXED-LAYER PARCEL, THEN
!                 LOOP OVER SUBSEQUENT DEPARTURE LAYERS IN LOWEST 350 hPa
!                 Theta_e =Theta*exp[L q_v/(C_p T)] 
!                         = T*(P0/P)**(R_d/C_p) * exp[L q_v/(C_p T)]
!                 -> THIS WILL BE THE UPDRAUGHT PARCEL (CONSERVING ITS
!                 PROPERTIES)  (no entrainment)
!                 CAPE    = Int ( g dTheta_v/Theta_v dz ) = 
!                   aprox = Int ( g (Theta_e_up-Theta_e_sat)/Theta_e_sat ) dz
!                 WITH THIS FORMULATION THE ACTUAL CAPE IS OVERESTIMATED  BY
!                 ROUGHLY 20%. DEEP CONVECTION CAN BE CONSIDERED FOR CAPE
!                 VALUES ABOVE 200-500 J/KG            
!
!     PARAMETER     DESCRIPTION                              UNITS
!     ---------     -----------                              -----
!     INPUT PARAMETERS (INTEGER):
!
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!
!    INPUT PARAMETERS (REAL):
!
!    *PAP*          PRESSURE ON FULL LEVELS                    PA
!    *PAPH*         PRESSURE ON HALF LEVELS                    PA
!    *PT*           TEMPERATURE ON FULL LEVELS                 K   
!    *PQ*           SPECIFIC HUMIDITY ON FULL LEVELS          KG/KG
!
!    OUTPUT PARAMETERS (REAL):
!
!    *PCAPE*        CONVECTIVE AVAILABLE POT. ENERGY          J/KG
!    *PCIN*         CONVECTIVE INHIBITION                     J/KG
!
!          MODIFICATIONS
!          -------------
!     24-06-2011 :  Add CIN    P. Bechtold
!
!-------------------------------------------------------------------------------
!
  subroutine cuancape2(kidia,kfdia,klon,klev,pap,paph,pt,pq,pcape,pcin)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    real(rk8) , dimension(klon,klev) , intent(in) :: pap
    real(rk8) , dimension(klon,klev+1) , intent(in) :: paph
    real(rk8) , dimension(klon,klev) , intent(in) :: pt
    real(rk8) , dimension(klon,klev) , intent(in) :: pq
    real(rk8) , dimension(klon) , intent(out) :: pcape
    real(rk8) , dimension(klon) , intent(out) :: pcin

    integer(ik4) :: jl , jk , jkk
    ! integer(ik4) , dimension(klon) :: jdep
    real(rk8) , dimension(klon) :: zpmix , ztmix , zthmix , zqmix , ztheteu
    real(rk8) , dimension(klon) :: zcin2
    real(rk8) , dimension(klon,klev) :: zcape , zcin
    real(rk8) :: zdp , zthetes , zqs , zdz , ztemp , zrpap
    real(rk8) , dimension(klon,klev) :: zthetad
    logical , dimension(klon,klev) :: llpap
    logical , dimension(klon) :: llz

    !----------------------------------------------------------------------

    do jl = kidia , kfdia
      pcape(jl) = d_zero
      pcin(jl) = d_zero
      zcin2(jl) = -1000.0D0
    end do
    do jk = klev-1 , njkt2,-1
      do jl = kidia , kfdia
        llpap(jl,jk) = ( pap(jl,jk)>80.D2 )
        if ( llpap(jl,jk) ) then
          zrpap = d_one/pap(jl,jk)
          zqs = foeewm(pt(jl,jk))*zrpap
          zqs = max(1.D-8,zqs)
          zqs = zqs/(d_one-retv*zqs) ! small correction
          zthetes = pt(jl,jk)*(stdp*zrpap)**rovcp * &
                    exp( foeldcp(pt(jl,jk))*zqs/pt(jl,jk) ) 
          zthetad(jl,jk) = d_one/zthetes
        end if
      end do
    end do
    do jkk = klev-1 , njkt1 , -1
      do jl = kidia , kfdia
        zcape(jl,jkk) = d_zero
        zcin(jl,jkk) = d_zero
        if ( paph(jl,klev+1)-paph(jl,jkk-1) < 60.D2 ) then
          ztmix(jl) = d_zero
          zthmix(jl) = d_zero
          zqmix(jl) = d_zero
          zpmix(jl) = d_zero
          do jk = jkk+1 , jkk-1 , -1
            if ( zpmix(jl) < 30.D2 ) then
              zdp = paph(jl,jk+1)-paph(jl,jk)
              zpmix(jl) = zpmix(jl)+zdp
              zthmix(jl) = zthmix(jl)+pt(jl,jk)*zdp*(stdp/pap(jl,jk))**rovcp
              zqmix(jl) = zqmix(jl)+pq(jl,jk)*zdp
            end if
          end do
          zdp = d_one/zpmix(jl)
          zqmix(jl) = zqmix(jl)*zdp
          zpmix(jl) = paph(jl,jkk+2) - 0.5D0*zpmix(jl)
          zthmix(jl) = zthmix(jl)*zdp
          ztmix(jl) = zthmix(jl)*(zpmix(jl)/stdp)**rovcp
        else
          zqmix(jl) = pq(jl,jkk)
          zpmix(jl) = pap(jl,jkk)
          ztmix(jl) = pt(jl,jkk)
          zthmix(jl) = pt(jl,jkk)*(stdp/zpmix(jl))**rovcp
        end if
        ztheteu(jl) = zthmix(jl)*exp(foeldcp(ztmix(jl))*zqmix(jl)/ztmix(jl))
        llz(jl) = (paph(jl,klev+1)-paph(jl,jkk)) < 350.D2
      end do
      do jk = jkk , njkt2 , -1
        do jl = kidia , kfdia
          if ( llpap(jl,jk) .and. llz(jl) ) then
            zrpap = d_one/pap(jl,jk)
            ztemp = ztheteu(jl)*zthetad(jl,jk)-d_one
            zdz = (paph(jl,jk+1)-paph(jl,jk))*zrpap*rgas*pt(jl,jk) * &
                  (d_one+retv*pq(jl,jk))
            if ( ztemp > d_zero ) then
              zcape(jl,jkk) = zcape(jl,jkk)+ztemp*zdz
            else if ( ztemp < d_zero .and. zcape(jl,jkk) < 100.D0 ) then
              zcin(jl,jkk) = zcin(jl,jkk)+ztemp*zdz
            end if
          end if
        end do
      end do
    end do
    ! chose maximum cape and cin values
    ! jdep(:) = klev-1
    do jk = klev-1 , njkt1 , -1
      do jl = kidia , kfdia
        if ( zcape(jl,jk) > pcape(jl) ) then
          pcape(jl) = zcape(jl,jk)
          ! jdep(jl) = jk
        end if
        if ( zcin(jl,jk) > zcin2(jl) .and. zcin(jl,jk) < d_zero ) then
          zcin2(jl) = zcin(jl,jk)
        end if
      end do
    end do
    do jl = kidia , kfdia
      ! jkk = jdep(jl)
      ! pcin(jl) = -zcin(jl,jkk)
      pcin(jl) = -zcin2(jl)
    end do
  end subroutine cuancape2
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
  subroutine cuentr(kidia,kfdia,klon,klev,kk,kcbot,ldcum,ldwork, &
                    pgeoh,pmfu,pdmfen,pdmfde)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    integer(ik4) , intent(in) :: kk
    integer(ik4) , dimension(klon) , intent(in) :: kcbot
    logical , dimension(klon) , intent(in) :: ldcum
    logical , intent(in) :: ldwork
    real(rk8) , dimension(klon,klev+1) , intent(in) :: pgeoh
    real(rk8) , dimension(klon,klev) , intent(in) :: pmfu
    real(rk8) , dimension(klon) , intent(out) :: pdmfen
    real(rk8) , dimension(klon) , intent(out) :: pdmfde
    logical :: llo1
    integer(ik4) :: jl
    real(rk8) :: zdz , zmf
    real(rk8) , dimension(klon) :: zentr
    !
    !* 1. CALCULATE ENTRAINMENT AND DETRAINMENT RATES
    ! -------------------------------------------
    if ( ldwork ) then
      do jl = kidia , kfdia
        pdmfen(jl) = d_zero
        pdmfde(jl) = d_zero
        zentr(jl) = d_zero
      end do
      !
      !*  1.1 SPECIFY ENTRAINMENT RATES
      !   -------------------------
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
  subroutine cuinin(kidia,kfdia,klon,klev,pten,pqen,pqsen,puen,pven,  &
                    pvervel,pgeo,paph,klwmin,klab,ptenh,pqenh,pqsenh, &
                    pgeoh,ptu,pqu,ptd,pqd,puu,pvu,pud,pvd,plu)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    real(rk8) , dimension(klon,klev) , intent(in) :: pten
    real(rk8) , dimension(klon,klev) , intent(in) :: pqen
    real(rk8) , dimension(klon,klev) , intent(in) :: pqsen
    real(rk8) , dimension(klon,klev) , intent(in) :: puen
    real(rk8) , dimension(klon,klev) , intent(in) :: pven
    real(rk8) , dimension(klon,klev) , intent(in) :: pvervel
    real(rk8) , dimension(klon,klev) , intent(in) :: pgeo
    real(rk8) , dimension(klon,klev+1) , intent(in) :: paph
    integer(ik4) , dimension(klon) , intent(out) :: klwmin
    integer(ik4) , dimension(klon,klev) , intent(out) :: klab
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptenh
    real(rk8) , dimension(klon,klev) , intent(out) :: pqenh
    real(rk8) , dimension(klon,klev) , intent(inout) :: pqsenh
    real(rk8) , dimension(klon,klev+1) , intent(in) :: pgeoh
    real(rk8) , dimension(klon,klev) , intent(out) :: ptu
    real(rk8) , dimension(klon,klev) , intent(out) :: pqu
    real(rk8) , dimension(klon,klev) , intent(out) :: ptd
    real(rk8) , dimension(klon,klev) , intent(out) :: pqd
    real(rk8) , dimension(klon,klev) , intent(out) :: puu
    real(rk8) , dimension(klon,klev) , intent(out) :: pvu
    real(rk8) , dimension(klon,klev) , intent(out) :: pud
    real(rk8) , dimension(klon,klev) , intent(out) :: pvd
    real(rk8) , dimension(klon,klev) , intent(out) :: plu
    real(rk8) , dimension(klon) :: zwmax
    real(rk8) , dimension(klon) :: zph
    logical , dimension(klon) :: llflag
    integer(ik4) :: icall , ik , jk , jl
    real(rk8) :: zalfa , zzs
    !----------------------------------------------------------------------
    ! 
    !*    1. SPECIFY LARGE SCALE PARAMETERS AT HALF LEVELS
    !*       ADJUST TEMPERATURE FIELDS IF STATICLY UNSTABLE
    !*       FIND LEVEL OF MAXIMUM VERTICAL VELOCITY
    ! ----------------------------------------------
    zalfa = log(d_two)
    do jk = 2 , klev
      do jl = kidia , kfdia
        ptenh(jl,jk) = (max(rcpd*pten(jl,jk-1) + &
          pgeo(jl,jk-1),rcpd*pten(jl,jk)+pgeo(jl,jk))-pgeoh(jl,jk))*rcpd
        pqenh(jl,jk) = pqen(jl,jk-1)
        pqsenh(jl,jk) = pqsen(jl,jk-1)
        zph(jl) = paph(jl,jk)
        llflag(jl) = .true.
      end do
      if ( jk >= klev-1 .or. jk < njkt2 ) cycle
      ik = jk
      icall = 3
      call cuadjtq(kidia,kfdia,klon,klev,ik,zph,ptenh,pqsenh,llflag,icall)
      do jl = kidia , kfdia
        pqenh(jl,jk) = min(pqen(jl,jk-1),pqsen(jl,jk-1)) + &
                          (pqsenh(jl,jk)-pqsen(jl,jk-1))
        pqenh(jl,jk) = max(pqenh(jl,jk),d_zero)
      end do
    end do
    do jl = kidia , kfdia
      ptenh(jl,klev) = (rcpd*pten(jl,klev)+pgeo(jl,klev)-pgeoh(jl,klev))*rcpd
      pqenh(jl,klev) = pqen(jl,klev)
      ptenh(jl,1) = pten(jl,1)
      pqenh(jl,1) = pqen(jl,1)
      klwmin(jl) = klev
      zwmax(jl) = d_zero
    end do
    do jk = klev - 1 , 2 , -1
      do jl = kidia , kfdia
        zzs = max(rcpd*ptenh(jl,jk)+pgeoh(jl,jk), &
                  rcpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))
        ptenh(jl,jk) = (zzs-pgeoh(jl,jk))*rcpd
      end do
    end do
    do jk = klev , 3 , -1
!DIR$ IVDEP
!OCL NOVREC
      do jl = kidia , kfdia
        if ( pvervel(jl,jk)<zwmax(jl) ) then
          zwmax(jl) = pvervel(jl,jk)
          klwmin(jl) = jk
        end if
      end do
    end do
    !-----------------------------------------------------------------------
    !*    2.0 INITIALIZE VALUES FOR UPDRAFTS AND DOWNDRAFTS
    !*        ---------------------------------------------
    do jk = 1 , klev
      ik = jk - 1
      if ( jk==1 ) ik = 1
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
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KLWMIN*       LEVEL OF MAXIMUM VERTICAL VELOCITY
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KCBOT*        CLOUD BASE LEVEL
!    *KDPL*         DEPARTURE LEVEL FOR CONVECTION
!
!    INPUT PARAMETERS (REAL):
!
!    *PTSPHY*       TIME STEP FOR THE PHYSICS                      S
!    *PTENH*        ENV. TEMPERATURE (T+1) ON HALF LEVELS          K
!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS     KG/KG
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)      M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)      M/S
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)      K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1) KG/KG
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)  KG/KG
!    *PGEO*         GEOPOTENTIAL                                 M2/S2
!    *PLITOT*       GRID MEAN LIQUID WATER+ICE CONTENT           KG/KG
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                  M2/S2
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS           PA
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
!    *PVERVEL*      VERTICAL VELOCITY                            PA/S
!
!    INPUT PARAMETERS (LOGICAL):
!
!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!
!    UPDATED PARAMETERS (INTEGER):
!
!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CLOUD LEVELS
!
!    UPDATED PARAMETERS (REAL):
!
!    *PTU*          TEMPERATURE IN UPDRAFTS                        K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                   KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS             KG/KG
!
!    OUTPUT PARAMETERS (INTEGER):
!
!    *KCTOP*        CLOUD TOP LEVEL
!    *KCTOP0*       FIRST GUESS OF CLOUD TOP LEVEL 
!
!    OUTPUT PARAMETERS (REAL):
!
!    *PMFU*         MASSFLUX IN UPDRAFTS                         KG/(M2*S)
!    *PMFUB*        MASSFLUX IN UPDRAFTS AT CLOUD BASE           KG/(M2*S)
!    *PMFUS*        FLUX OF DRY STATIC ENERGY IN UPDRAFTS         J/(M2*S)
!    *PMFUQ*        FLUX OF SPEC. HUMIDITY IN UPDRAFTS           KG/(M2*S)
!    *PMFUL*        FLUX OF LIQUID WATER IN UPDRAFTS             KG/(M2*S)
!    *PLUDE*        DETRAINED LIQUID WATER                       KG/(M2*S)
!    *PLGLAC*       FROZEN CLOUD WATER CONTENT                   KG/KG
!    *PDMFUP*       FLUX DIFFERENCE OF PRECIP. IN UPDRAFTS       KG/(M2*S)
!    *PMFUDE_RATE*  UPDRAFT DETRAINMENT RATE                     KG/(M2*S)
!    *PKINEU*       UPDRAFT KINETIC ENERGY                       M2/S2
!    *PWMEAN*       MEAN UPDRAUGHT VELOCITY                      M/S
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
!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             99-06-14 : Optimisation        D.SALMOND
!             01-05-22 : Modified flux limiter M.CULLEN
!             02-08-14 : Allow for departure level =/ KLEV  P.BECHTOLD
!             03-08-28 : Clean-up detrainment rates         P.BECHTOLD
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        J.Hague       08-Dec-2005 Tuning: LLFLAG indexing
!             07-06-01 : Organized entrainment based on RH  P.BECHTOLD
!
!----------------------------------------------------------------------
!
  subroutine cuascn(kidia,kfdia,klon,klev,ptsphy,ptenh,pqenh,pten,pqen, &
                    pqsen,plitot,pgeo,pgeoh,pap,paph,pvervel,pwubase,   &
                    ldland,ldcum,ktype,klab,ptu,pqu,plu,pmfu,pmfub,     &
                    plglac,pmfus,pmfuq,pmful,plude,pdmfup,pdmfen,kcbot, &
                    kctop,kctop0,kdpl,pmfude_rate,pkineu,pwmean)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    real(rk8) , intent(in) :: ptsphy
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptenh
    real(rk8) , dimension(klon,klev) , intent(inout) :: pqenh
    real(rk8) , dimension(klon,klev) , intent(in) :: pten
    real(rk8) , dimension(klon,klev) , intent(in) :: pqen
    real(rk8) , dimension(klon,klev) , intent(in) :: pqsen
    real(rk8) , dimension(klon,klev) , intent(in) :: plitot
    real(rk8) , dimension(klon,klev) , intent(in) :: pgeo
    real(rk8) , dimension(klon,klev+1) , intent(in) :: pgeoh
    real(rk8) , dimension(klon,klev) , intent(in) :: pap
    real(rk8) , dimension(klon,klev+1) , intent(in) :: paph
    real(rk8) , dimension(klon,klev) , intent(in) :: pvervel
    real(rk8) , dimension(klon) , intent(in) :: pwubase
    logical , dimension(klon) , intent(in) :: ldland
    logical , dimension(klon) , intent(inout) :: ldcum
    integer(ik4) , dimension(klon) , intent(inout) :: ktype
    integer(ik4) , dimension(klon,klev) , intent(inout) :: klab
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptu
    real(rk8) , dimension(klon,klev) , intent(inout) :: pqu
    real(rk8) , dimension(klon,klev) , intent(inout) :: plu
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfu
    real(rk8) , dimension(klon) , intent(inout) :: pmfub
    real(rk8) , dimension(klon,klev) , intent(out) :: plglac
    real(rk8) , dimension(klon,klev) , intent(out) :: pmfus
    real(rk8) , dimension(klon,klev) , intent(out) :: pmfuq
    real(rk8) , dimension(klon,klev) , intent(out) :: pmful
    real(rk8) , dimension(klon,klev) , intent(out) :: plude
    real(rk8) , dimension(klon,klev) , intent(out) :: pdmfup
    real(rk8) , dimension(klon,klev) , intent(out) :: pdmfen
    integer(ik4) , dimension(klon) , intent(inout) :: kcbot
    integer(ik4) , dimension(klon) , intent(out) :: kctop
    integer(ik4) , dimension(klon) , intent(inout) :: kctop0
    integer(ik4) , dimension(klon) , intent(in) :: kdpl
    real(rk8) , dimension(klon,klev) , intent(out) :: pmfude_rate
    real(rk8) , dimension(klon,klev) , intent(out) :: pkineu
    real(rk8) , dimension(klon) , intent(out) :: pwmean
    real(rk8) , dimension(klon,klev) :: zlrain , zbuo
    real(rk8) , dimension(klon) :: zdmfen , zdmfde , zqold , zluold , zprecip
    real(rk8) , dimension(klon) :: zdpmean
    real(rk8) , dimension(klon) :: zoentr , zph
    logical , dimension(klon) :: llflag , llflaguv , llo1
    logical :: llo3 , llo4
    integer(ik4) :: icall , ik , is , jk , jl , ikb
    integer(ik4) :: jll , jlm
    integer(ik4) , dimension(klon) :: jlx
    real(rk8) :: z_cldmax , z_cprc2 , z_cwdrag , z_cwifrac , zalfaw , zbc ,&
                 zbe , zbuoc , zc , zcbf , zcons2 , zd , zdfi , zdkbuo ,   &
                 zdken , zdnoprc , zdt , zfac , zfacbuo , zint ,           &
                 zkedke , zlcrit , zleen , zlnew , zmfmax , zmftest ,      &
                 zmfulk , zmfun , zmfuqk , zmfusk , zprcdgw , zprcon ,     &
                 zqeen , zqude , zrnew , zrold , zscde , zseen , ztglace , &
                 zvi , zvv , zvw , zwu , zzco
    real(rk8) :: zchange , zxs , zxe
    logical , dimension(klon) :: llklab
    !----------------------------------------------------------------------
    !*    1.           SPECIFY PARAMETERS
    ! ------------------
    zcons2 = rmfcfl/(egrav*ptsphy)
    ztglace = tzero - 13.0D0
    zfacbuo = d_half/(d_one+d_half)
    zprcdgw = rprcon/egrav
    z_cldmax = 5.D-3
    z_cwifrac = d_half
    z_cprc2 = d_half
    z_cwdrag = (3.0D0/8.0D0)*0.506D0/0.200D0
    !----------------------------------------------------------------------
    ! 2.           SET DEFAULT VALUES
    ! ------------------
    llo3 = .false.
    do jl = kidia , kfdia
      zluold(jl) = d_zero
      if ( .not. ldcum(jl) ) then
        kcbot(jl) = -1
        pmfub(jl) = d_zero
        pqu(jl,klev) = d_zero
        ktype(jl) = 0
      end if
      pwmean(jl) = d_zero
      zdpmean(jl) = d_zero
      zoentr(jl) = d_zero
    end do
    ! initalize various quantities
    ! note that liquid water and kinetic energy at cloud base is
    ! preserved from cubase
    do jl = kidia , kfdia
      llklab(jl) = .false.
      if ( .not. ldcum(jl) .or. ktype(jl) == 3 ) llklab(jl) = .true.
    end do
    do jk = 1 , klev
      do jl = kidia , kfdia
        if ( jk /= kcbot(jl) ) plu(jl,jk) = d_zero
        pkineu(jl,jk) = d_zero
      end do
      do jl = kidia , kfdia
        pmfu(jl,jk) = d_zero
        pmfus(jl,jk) = d_zero
        pmfuq(jl,jk) = d_zero
        pmful(jl,jk) = d_zero
      end do
      do jl = kidia , kfdia
        plude(jl,jk) = d_zero
        plglac(jl,jk) = d_zero
        pdmfup(jl,jk) = d_zero
        zlrain(jl,jk) = d_zero
      end do
      do jl = kidia , kfdia
        zbuo(jl,jk) = d_zero
        if ( llklab(jl) ) klab(jl,jk) = 0
        if ( .not. ldcum(jl) .and. paph(jl,jk) < 4.0D4 ) kctop0(jl) = jk
        pdmfen(jl,jk) = d_zero
        pmfude_rate(jl,jk) = d_zero
      end do
    end do
!DIR$ IVDEP
!OCL NOVREC
    do jl = kidia , kfdia
      if ( ktype(jl) == 3 ) ldcum(jl) = .false.
    end do
    !----------------------------------------------------------------------
    ! 3.0          INITIALIZE VALUES AT cloud base LEVEL
    ! -------------------------------------
    do jl = kidia , kfdia
      kctop(jl) = kcbot(jl)
      if ( ldcum(jl) ) then
        ikb = kcbot(jl)
        pkineu(jl,ikb) = d_half*pwubase(jl)**2
        pmfu(jl,ikb) = pmfub(jl)
        pmfus(jl,ikb) = pmfub(jl)*(rcpd*ptu(jl,ikb)+pgeoh(jl,ikb))
        pmfuq(jl,ikb) = pmfub(jl)*pqu(jl,ikb)
        pmful(jl,ikb) = pmfub(jl)*plu(jl,ikb)
      end if
    end do
    !----------------------------------------------------------------------
    ! 4.           DO ASCENT: SUBCLOUD LAYER (KLAB=1) ,CLOUDS (KLAB=2)
    ! BY DOING FIRST DRY-ADIABATIC ASCENT AND THEN
    ! BY ADJUSTING T,Q AND L ACCORDINGLY IN *CUADJTQ*,
    ! THEN CHECK FOR BUOYANCY AND SET FLAGS ACCORDINGLY
    ! -------------------------------------------------
    do jk = klev - 1 , 3 , -1
      !   SPECIFY CLOUD BASE VALUES FOR MIDLEVEL CONVECTION
      !   IN *CUBASMC* IN CASE THERE IS NOT ALREADY CONVECTION
      !   ----------------------------------------------------
      ik = jk
      call cubasmcn(kidia,kfdia,klon,klev,ik,pten,pqen,pqsen,  &
                    pvervel,pgeo,pgeoh,ldcum,ktype,klab,kcbot, &
                    pmfu,pmfub,zlrain,ptu,pqu,plu,pmfus,pmfuq, &
                    pmful,pdmfup)
      is = 0
      jlm = 0
      do jl = kidia , kfdia
        llflag(jl) = .false.
        zprecip(jl) = d_zero
        llo1(jl) = .false.
        is = is + klab(jl,jk+1)
        if ( klab(jl,jk+1) == 0 ) klab(jl,jk) = 0
        if ( (ldcum(jl) .and. klab(jl,jk+1) == 2) .or.   &
             (ktype(jl) == 3 .and. klab(jl,jk+1) == 1) ) then
          llflag(jl) = .true.
          jlm = jlm + 1
          jlx(jlm) = jl
        end if
        if ( klab(jl,jk+1) > 0 ) then
          llflaguv(jl) = .true.
        else
          llflaguv(jl) = .false.
        end if
        zph(jl) = paph(jl,jk)
        if ( ktype(jl) == 3 .and. jk == kcbot(jl) ) then
          zmfmax = (paph(jl,jk)-paph(jl,jk-1))*zcons2*rmflic + rmflia
          if ( pmfub(jl) > zmfmax ) then
            zfac = zmfmax/pmfub(jl)
            pmfu(jl,jk+1) = pmfu(jl,jk+1)*zfac
            pmfus(jl,jk+1) = pmfus(jl,jk+1)*zfac
            pmfuq(jl,jk+1) = pmfuq(jl,jk+1)*zfac
            pmfub(jl) = zmfmax
          end if
        end if
      end do
      if ( is > 0 ) llo3 = .true.
      !*  SPECIFY ENTRAINMENT RATES IN *CUENTR*
      !   -------------------------------------
      ik = jk
      call cuentr(kidia,kfdia,klon,klev,ik,kcbot,ldcum,llo3, &
                  pgeoh,pmfu,zdmfen,zdmfde)
      !   DO ADIABATIC ASCENT FOR ENTRAINING/DETRAINING PLUME
      !   ---------------------------------------------------
      if ( llo3 ) then
        llo4 = ptsphy > 1800.0D0 .and. rmfcfl == d_one
        do jl = kidia , kfdia
          zqold(jl) = d_zero
        end do
        do jll = 1 , jlm
          jl = jlx(jll)
          zdmfde(jl) = min(zdmfde(jl),0.75D0*pmfu(jl,jk+1))
          if ( jk == kcbot(jl) ) then
            zoentr(jl) = -entrorg*(min(d_one,pqen(jl,jk)/pqsen(jl,jk)) - &
                         d_one)*(pgeoh(jl,jk)-pgeoh(jl,jk+1))*regrav
            zoentr(jl) = min(0.4D0,zoentr(jl))*pmfu(jl,jk+1)
          end if
          if ( jk < kcbot(jl) ) then
            zmfmax = (paph(jl,jk)-paph(jl,jk-1))*zcons2*rmflic + rmflia
            if ( ktype(jl) == 2 .and. llo4 ) zmfmax = zmfmax*3.0D0
            zxs = max(pmfu(jl,jk+1)-zmfmax,d_zero)
            pwmean(jl) = pwmean(jl) + pkineu(jl,jk+1)*(pap(jl,jk+1)-pap(jl,jk))
            zdpmean(jl) = zdpmean(jl) + pap(jl,jk+1) - pap(jl,jk)
            zdmfen(jl) = zoentr(jl)
            if ( ktype(jl) >= 2 ) then
              zdmfen(jl) = entshalp*zdmfen(jl)
              zdmfde(jl) = zdmfen(jl)
            end if
            zdmfde(jl) = zdmfde(jl) * &
                         (1.6D0-min(d_one,pqen(jl,jk)/pqsen(jl,jk)))
            zmftest = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
            zchange = max(zmftest-zmfmax,d_zero)
            zxe = max(zchange-zxs,d_zero)
            zdmfen(jl) = zdmfen(jl) - zxe
            zchange = zchange - zxe
            zdmfde(jl) = zdmfde(jl) + zchange
          end if
          pdmfen(jl,jk) = zdmfen(jl) - zdmfde(jl)
          pmfu(jl,jk) = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
          zqeen = pqenh(jl,jk+1)*zdmfen(jl)
          zseen = (rcpd*ptenh(jl,jk+1)+pgeoh(jl,jk+1))*zdmfen(jl)
          if ( plitot(jl,jk)>rlmin ) then
            zleen = plitot(jl,jk)*zdmfen(jl)
          else
            zleen = d_zero
          end if
          zscde = (rcpd*ptu(jl,jk+1)+pgeoh(jl,jk+1))*zdmfde(jl)
          zqude = pqu(jl,jk+1)*zdmfde(jl)
          plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
          zmfusk = pmfus(jl,jk+1) + zseen - zscde
          zmfuqk = pmfuq(jl,jk+1) + zqeen - zqude
          zmfulk = pmful(jl,jk+1) + zleen - plude(jl,jk)
          plu(jl,jk) = zmfulk*(d_one/max(rmfcmin,pmfu(jl,jk)))
          pqu(jl,jk) = zmfuqk*(d_one/max(rmfcmin,pmfu(jl,jk)))
          ptu(jl,jk) = (zmfusk * &
            (d_one/max(rmfcmin,pmfu(jl,jk)))-pgeoh(jl,jk))/rcpd
          ptu(jl,jk) = max(100.0D0,ptu(jl,jk))
          ptu(jl,jk) = min(400.0D0,ptu(jl,jk))
          zqold(jl) = pqu(jl,jk)
          zlrain(jl,jk) = zlrain(jl,jk+1)*(pmfu(jl,jk+1)-zdmfde(jl)) * &
                          (d_one/max(rmfcmin,pmfu(jl,jk)))
          zluold(jl) = plu(jl,jk)
        end do
        ! reset to environmental values if below departure level
        do jl = kidia , kfdia
          if ( jk > kdpl(jl) ) then
            ptu(jl,jk) = ptenh(jl,jk)
            pqu(jl,jk) = pqenh(jl,jk)
            plu(jl,jk) = d_zero
            zluold(jl) = plu(jl,jk)
          end if
        end do
        ! DO CORRECTIONS FOR MOIST ASCENT
        ! BY ADJUSTING T,Q AND L IN *CUADJTQ*
        ! -----------------------------------
        ik = jk
        icall = 1
        if ( jlm > 0 ) then
          call cuadjtq(kidia,kfdia,klon,klev,ik,zph,ptu,pqu,llflag,icall)
        end if
!DIR$   IVDEP
!OCL    NOVREC
        do jll = 1 , jlm
          jl = jlx(jll)
          if ( pqu(jl,jk) /= zqold(jl) ) then
            plglac(jl,jk) = plu(jl,jk) * &
                           ((d_one-foealfcu(ptu(jl,jk)))- &
                            (d_one-foealfcu(ptu(jl,jk+1))))
            ptu(jl,jk) = ptu(jl,jk) + wlhfocp*plglac(jl,jk)
          end if
        end do
        do jll = 1 , jlm
          jl = jlx(jll)
          if ( pqu(jl,jk) /= zqold(jl) ) then
            klab(jl,jk) = 2
            plu(jl,jk) = plu(jl,jk) + zqold(jl) - pqu(jl,jk)
            zbc = ptu(jl,jk)*(d_one+retv*pqu(jl,jk)-plu(jl,jk+1) - &
              zlrain(jl,jk+1))
            zbe = ptenh(jl,jk)*(d_one+retv*pqenh(jl,jk))
            zbuo(jl,jk) = zbc - zbe
            ! set flags in case of midlevel convection
            if ( ktype(jl) == 3 .and. klab(jl,jk+1) == 1 ) then
              if ( zbuo(jl,jk) > -d_half ) then
                ldcum(jl) = .true.
                kctop(jl) = jk
                pkineu(jl,jk) = d_half
              else
                klab(jl,jk) = 0
                pmfu(jl,jk) = d_zero
                plude(jl,jk) = d_zero
                plu(jl,jk) = d_zero
              end if
            end if
            if ( klab(jl,jk+1) == 2 ) then
              if ( zbuo(jl,jk) < d_zero .and. klab(jl,jk+1) == 2 ) then
                ptenh(jl,jk) = d_half*(pten(jl,jk)+pten(jl,jk-1))
                pqenh(jl,jk) = d_half*(pqen(jl,jk)+pqen(jl,jk-1))
                zbuo(jl,jk) = zbc - ptenh(jl,jk)*(d_one+retv*pqenh(jl,jk))
              end if
              zbuoc = (zbuo(jl,jk) / &
                (ptenh(jl,jk)*(d_one+retv*pqenh(jl,jk)))+zbuo(jl,jk+1) / &
                (ptenh(jl,jk+1)*(d_one+retv*pqenh(jl,jk+1))))*d_half
              zdkbuo = (pgeoh(jl,jk)-pgeoh(jl,jk+1))*zfacbuo*zbuoc
              ! mixing and "pressure" gradient term in upper
              ! troposphere
              if ( zdmfen(jl) > d_zero ) then
                zdken = min(d_one,(d_one+z_cwdrag)*zdmfen(jl) / &
                        max(rmfcmin,pmfu(jl,jk+1)))
              else
                zdken = min(d_one,(d_one+z_cwdrag)*zdmfde(jl) / &
                        max(rmfcmin,pmfu(jl,jk+1)))
              end if
              pkineu(jl,jk) = (pkineu(jl,jk+1)*(d_one-zdken)+zdkbuo) / &
                (d_one+zdken)
              if ( zbuo(jl,jk) < d_zero .and. klab(jl,jk+1) == 2 ) then
                zkedke = pkineu(jl,jk)/max(1.D-10,pkineu(jl,jk+1))
                zkedke = max(d_zero,min(d_one,zkedke))
                zmfun = sqrt(zkedke)*pmfu(jl,jk+1)
                zdmfde(jl) = max(zdmfde(jl),pmfu(jl,jk+1)-zmfun)
                plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
                pmfu(jl,jk) = pmfu(jl,jk+1) + zdmfen(jl) - zdmfde(jl)
              end if
              if ( zbuo(jl,jk) >- 0.2D0 .and. klab(jl,jk+1) == 2 ) then
                ikb = kcbot(jl)
                zoentr(jl) = entrorg*(0.3D0-(min(d_one,pqen(jl,jk-1) /    &
                  pqsen(jl,jk-1))-d_one))*(pgeoh(jl,jk-1)-pgeoh(jl,jk)) * &
                  regrav*min(d_one,pqsen(jl,jk)/pqsen(jl,ikb))**3
                zoentr(jl) = min(0.4D0,zoentr(jl))*pmfu(jl,jk)
              else
                zoentr(jl) = d_zero
              end if
             ! Erase values if below departure level
              if ( jk > kdpl(jl) ) then
                pmfu(jl,jk) = pmfu(jl,jk+1)
                pkineu(jl,jk) = d_half
              end if
              if ( pkineu(jl,jk) > d_zero .and. pmfu(jl,jk) > d_zero ) then
                kctop(jl) = jk
                llo1(jl) = .true.
              else
                klab(jl,jk) = 0
                pmfu(jl,jk) = d_zero
                pkineu(jl,jk) = d_zero
                zdmfde(jl) = pmfu(jl,jk+1)
                plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
              end if
              ! store detrainment rates for updraught
              if ( pmfu(jl,jk+1) > d_zero ) pmfude_rate(jl,jk) = zdmfde(jl)
            end if
          else if ( ktype(jl) == 2 .and. pqu(jl,jk) == zqold(jl) ) then
            klab(jl,jk) = 0
            pmfu(jl,jk) = d_zero
            pkineu(jl,jk) = d_zero
            zdmfde(jl) = pmfu(jl,jk+1)
            plude(jl,jk) = plu(jl,jk+1)*zdmfde(jl)
            pmfude_rate(jl,jk) = zdmfde(jl)
          end if
        end do
        !     CALCULATE PRECIPITATION RATE BY
        !     ANALYTIC INTEGRATION OF EQUATION FOR L
        do jl = kidia , kfdia
          if ( llo1(jl) ) then
            if ( ldland(jl) ) then
              zdnoprc = 5.D-4
            else
              zdnoprc = 3.D-4
            end if
            if ( plu(jl,jk) > zdnoprc ) then
              zwu = min(15.0D0,sqrt(d_two*max(0.1D0,pkineu(jl,jk+1))))
              zprcon = zprcdgw/(0.75D0*zwu)
              ! PARAMETERS FOR BERGERON-FINDEISEN PROCESS (T < -5C)
              zdt = min(rtber-rtice,max(rtber-ptu(jl,jk),d_zero))
              zcbf = d_one + z_cprc2*sqrt(zdt)
              zzco = zprcon*zcbf
              zlcrit = zdnoprc/zcbf
              zdfi = pgeoh(jl,jk) - pgeoh(jl,jk+1)
              zc = (plu(jl,jk)-zluold(jl))
              zd = zzco*(d_one-exp(-(plu(jl,jk)/zlcrit)**2))*zdfi
              zint = exp(-zd)
              zlnew = zluold(jl)*zint + zc/zd*(d_one-zint)
              zlnew = max(d_zero,min(plu(jl,jk),zlnew))
              zlnew = min(z_cldmax,zlnew)
              zprecip(jl) = max(d_zero,zluold(jl)+zc-zlnew)
              pdmfup(jl,jk) = zprecip(jl)*pmfu(jl,jk)
              zlrain(jl,jk) = zlrain(jl,jk) + zprecip(jl)
              plu(jl,jk) = zlnew
            end if
          end if
        end do
        do jl = kidia , kfdia
          if ( llo1(jl) ) then
            if ( zlrain(jl,jk) > d_zero ) then
              zvw = 21.18D0*zlrain(jl,jk)**0.2D0
              zvi = z_cwifrac*zvw
              zalfaw = foealfcu(ptu(jl,jk))
              zvv = zalfaw*zvw + (d_one-zalfaw)*zvi
              zrold = zlrain(jl,jk) - zprecip(jl)
              zc = zprecip(jl)
              zwu = min(15.0D0,sqrt(d_two*max(0.1D0,pkineu(jl,jk))))
              zd = zvv/zwu
              zint = exp(-zd)
              zrnew = zrold*zint + zc/zd*(d_one-zint)
              zrnew = max(d_zero,min(zlrain(jl,jk),zrnew))
              zlrain(jl,jk) = zrnew
            end if
          end if
        end do
        do jll = 1 , jlm
          jl = jlx(jll)
          pmful(jl,jk) = plu(jl,jk)*pmfu(jl,jk)
          pmfus(jl,jk) = (rcpd*ptu(jl,jk)+pgeoh(jl,jk))*pmfu(jl,jk)
          pmfuq(jl,jk) = pqu(jl,jk)*pmfu(jl,jk)
        end do
      end if
    end do
    !----------------------------------------------------------------------
    ! 5.           FINAL CALCULATIONS
    ! ------------------
    do jl = kidia , kfdia
      if ( kctop(jl) == -1 ) ldcum(jl) = .false.
      kcbot(jl) = max(kcbot(jl),kctop(jl))
      if ( ldcum(jl) ) then
        pwmean(jl) = max(1.D-2,pwmean(jl)/max(d_one,zdpmean(jl)))
        pwmean(jl) = sqrt(d_two*pwmean(jl))
      end if
    end do
  end subroutine cuascn
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
  subroutine cuadjtq(kidia,kfdia,klon,klev,kk,psp,pt,pq,ldflag,kcall)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    integer(ik4) , intent(in) :: kk
    real(rk8) , dimension(klon) , intent(in) :: psp
    real(rk8) , dimension(klon,klev) , intent(inout) :: pt
    real(rk8) , dimension(klon,klev) , intent(inout) :: pq
    logical , dimension(klon) , intent(in) :: ldflag
    integer(ik4) , intent(in) :: kcall
    integer(ik4) :: jl , jlen

    real(rk8) , dimension(kfdia-kidia+1) :: ztmp0
    real(rk8) , dimension(kfdia-kidia+1) :: ztmp1
    real(rk8) , dimension(kfdia-kidia+1) :: ztmp2
    real(rk8) , dimension(kfdia-kidia+1) :: ztmp3
    real(rk8) , dimension(kfdia-kidia+1) :: ztmp4
    real(rk8) , dimension(kfdia-kidia+1) :: ztmp5
    real(rk8) , dimension(kfdia-kidia+1) :: ztmp6

    real(rk8) :: zcond , zcond1 , zcor , zqmax , zqsat , zqp
    real(rk8) :: zl , zi , zf

    !----------------------------------------------------------------------
    ! 1.           DEFINE CONSTANTS
    ! ----------------

    if ( n_vmass > 0 ) jlen = kfdia - kidia + 1

    zqmax = d_half

    !   2.           CALCULATE CONDENSATION AND ADJUST T AND Q ACCORDINGLY
    !   -----------------------------------------------------


    if ( kcall == 1 ) then
!DIR$ IVDEP
!OCL  NOVREC
      do jl = kidia , kfdia
        if ( ldflag(jl) ) then
          zqp = d_one/psp(jl)
          zl = d_one/(pt(jl,kk)-c4les)
          zi = d_one/(pt(jl,kk)-c4ies)
          zqsat = c2es*(foealfcu(pt(jl,kk))*exp(c3les*(pt(jl,kk)-tzero)*zl) + &
                (d_one-foealfcu(pt(jl,kk)))*exp(c3ies*(pt(jl,kk)-tzero)*zi))
          zqsat = zqsat*zqp
          zqsat = min(d_half,zqsat)
          zcor = d_one - retv*zqsat
          zf = foealfcu(pt(jl,kk))*c5alvcp*zl**2 + &
               (d_one-foealfcu(pt(jl,kk)))*c5alscp*zi**2
          zcond = (pq(jl,kk)*zcor**2-zqsat*zcor)/(zcor**2+zqsat*zf)
          if ( zcond > d_zero ) then
            pt(jl,kk) = pt(jl,kk) + foeldcpmcu(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk) - zcond
            zl = d_one/(pt(jl,kk)-c4les)
            zi = d_one/(pt(jl,kk)-c4ies)
            zqsat = c2es*(foealfcu(pt(jl,kk)) * &
              exp(c3les*(pt(jl,kk)-tzero)*zl)+(d_one-foealfcu(pt(jl,kk))) * &
              exp(c3ies*(pt(jl,kk)-tzero)*zi))
            zqsat = zqsat*zqp
            zqsat = minj(d_half,zqsat)
            zcor = d_one - retv*zqsat
            zf = foealfcu(pt(jl,kk))*c5alvcp*zl**2 + &
                 (d_one-foealfcu(pt(jl,kk)))*c5alscp*zi**2
            zcond1 = (pq(jl,kk)*zcor**2-zqsat*zcor)/(zcor**2+zqsat*zf)
            if ( zcond == d_zero ) zcond1 = d_zero
            pt(jl,kk) = pt(jl,kk) + foeldcpmcu(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk) - zcond1
          end if
        end if
      end do
    else if ( kcall == 2 ) then
!DIR$ IVDEP
!OCL  NOVREC
      do jl = kidia , kfdia
        if ( ldflag(jl) ) then
          zqp = d_one/psp(jl)
          zqsat = foeewmcu(pt(jl,kk))*zqp
          zqsat = min(d_half,zqsat)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedemcu(pt(jl,kk)))
          zcond = min(zcond,d_zero)
          pt(jl,kk) = pt(jl,kk) + foeldcpmcu(pt(jl,kk))*zcond
          pq(jl,kk) = pq(jl,kk) - zcond
          zqsat = foeewmcu(pt(jl,kk))*zqp
          zqsat = min(d_half,zqsat)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedemcu(pt(jl,kk)))
          if ( zcond == d_zero ) zcond1 = min(zcond1,d_zero)
          pt(jl,kk) = pt(jl,kk) + foeldcpmcu(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk) - zcond1
        end if
      end do
    else if ( kcall == 0 ) then
!DIR$ IVDEP
!OCL  NOVREC
      do jl = kidia , kfdia
        zqp = d_one/psp(jl)
        zqsat = foeewm(pt(jl,kk))*zqp
        zqsat = min(d_half,zqsat)
        zcor = d_one/(d_one-retv*zqsat)
        zqsat = zqsat*zcor
        zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
        pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
        pq(jl,kk) = pq(jl,kk) - zcond1
        zqsat = foeewm(pt(jl,kk))*zqp
        zqsat = min(d_half,zqsat)
        zcor = d_one/(d_one-retv*zqsat)
        zqsat = zqsat*zcor
        zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
        pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
        pq(jl,kk) = pq(jl,kk) - zcond1
      end do
    else if ( kcall == 4 ) then
!DIR$ IVDEP
!OCL  NOVREC
      do jl = kidia , kfdia
        if ( ldflag(jl) ) then
          zqp = d_one/psp(jl)
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(d_half,zqsat)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond
          pq(jl,kk) = pq(jl,kk) - zcond
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(d_half,zqsat)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk) - zcond1
        end if
      end do
    else if ( kcall == 5 ) then ! Same as 4 but with LDFLAG all true
!DIR$ IVDEP
!OCL  NOVREC
      if ( n_vmass <= 0 ) then ! Not using Vector MASS
        do jl = kidia , kfdia
          zqp = d_one/psp(jl)
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(d_half,zqsat)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          zcond = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond
          pq(jl,kk) = pq(jl,kk) - zcond
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(d_half,zqsat)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedem(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk) - zcond1
        end do
      else ! Using Vector VMASS
        do jl = kidia , kfdia
          ztmp1(jl-kidia+1) = c3les*(pt(jl,kk)-tzero)
          ztmp2(jl-kidia+1) = c3ies*(pt(jl,kk)-tzero)
          ztmp3(jl-kidia+1) = pt(jl,kk) - c4les
          ztmp4(jl-kidia+1) = pt(jl,kk) - c4ies
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
          zqsat = minj(d_half,zqsat)
          zcor = d_one - retv*zqsat
          zf = foealfa(pt(jl,kk))*c5alvcp*(ztmp5(jl-kidia+1)**2) + &
               (d_one-foealfa(pt(jl,kk)))*c5alscp*(ztmp6(jl-kidia+1)**2)
          zcond = (pq(jl,kk)*zcor**2-zqsat*zcor)/(zcor**2+zqsat*zf)
          pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond
          pq(jl,kk) = pq(jl,kk) - zcond
          ztmp0(jl-kidia+1) = zqp
          ztmp1(jl-kidia+1) = c3les*(pt(jl,kk)-tzero)
          ztmp2(jl-kidia+1) = c3ies*(pt(jl,kk)-tzero)
          ztmp3(jl-kidia+1) = pt(jl,kk) - c4les
          ztmp4(jl-kidia+1) = pt(jl,kk) - c4ies
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
          zqsat = minj(d_half,zqsat)
          zcor = d_one - retv*zqsat
          zf = foealfa(pt(jl,kk))*c5alvcp*(ztmp5(jl-kidia+1)**2) + &
               (d_one-foealfa(pt(jl,kk)))*c5alscp*(ztmp6(jl-kidia+1)**2)
          zcond1 = (pq(jl,kk)*zcor**2-zqsat*zcor)/(zcor**2+zqsat*zf)
          pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk) - zcond1
        end do
      end if
    else if ( kcall == 3 ) then
!DIR$ IVDEP !OCL NOVREC
      if ( n_vmass <= 0 ) then ! Not using Vector MASS
        do jl = kidia , kfdia
          zqp = d_one/psp(jl)
          zqsat = foeewmcu(pt(jl,kk))*zqp
          zqsat = min(d_half,zqsat)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedemcu(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk) + foeldcpmcu(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk) - zcond1
          zqsat = foeewmcu(pt(jl,kk))*zqp
          zqsat = min(d_half,zqsat)
          zcor = d_one/(d_one-retv*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(d_one+zqsat*zcor*foedemcu(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk) + foeldcpmcu(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk) - zcond1
        end do
      else
        do jl = kidia , kfdia
          ztmp1(jl-kidia+1) = c3les*(pt(jl,kk)-tzero)
          ztmp2(jl-kidia+1) = c3ies*(pt(jl,kk)-tzero)
          ztmp3(jl-kidia+1) = pt(jl,kk) - c4les
          ztmp4(jl-kidia+1) = pt(jl,kk) - c4ies
        end do
        call vdiv(ztmp5,ztmp1,ztmp3,jlen)
        call vdiv(ztmp6,ztmp2,ztmp4,jlen)
        call vexp(ztmp1,ztmp5,jlen)
        call vexp(ztmp2,ztmp6,jlen)
        do jl = kidia , kfdia
          zqp = d_one/psp(jl)
          zqsat = c2es*(foealfcu(pt(jl,kk))*ztmp1(jl-kidia+1) + &
                 (d_one-foealfcu(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
          zqsat = minj(d_half,zqsat)
          zcor = d_one - retv*zqsat
          zcond1 = (pq(jl,kk)*zcor**2-zqsat*zcor) / &
                   (zcor**2+zqsat*foedemcu(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk) + foeldcpmcu(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk) - zcond1
          ztmp0(jl-kidia+1) = zqp
          ztmp1(jl-kidia+1) = c3les*(pt(jl,kk)-tzero)
          ztmp2(jl-kidia+1) = c3ies*(pt(jl,kk)-tzero)
          ztmp3(jl-kidia+1) = pt(jl,kk) - c4les
          ztmp4(jl-kidia+1) = pt(jl,kk) - c4ies
        end do
        call vdiv(ztmp5,ztmp1,ztmp3,jlen)
        call vdiv(ztmp6,ztmp2,ztmp4,jlen)
        call vexp(ztmp1,ztmp5,jlen)
        call vexp(ztmp2,ztmp6,jlen)
        do jl = kidia , kfdia
          zqp = ztmp0(jl-kidia+1)
          zqsat = c2es*(foealfcu(pt(jl,kk))*ztmp1(jl-kidia+1) + &
                  (d_one-foealfcu(pt(jl,kk)))*ztmp2(jl-kidia+1))*zqp
          zqsat = minj(d_half,zqsat)
          zcor = d_one - retv*zqsat
          zcond1 = (pq(jl,kk)*zcor**2-zqsat*zcor) / &
                   (zcor**2+zqsat*foedemcu(pt(jl,kk)))
          pt(jl,kk) = pt(jl,kk) + foeldcpmcu(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk) - zcond1
        end do
      end if
    else
      call fatal(__FILE__,__LINE__,'Unknown method kcall in cuadjtq')
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
  subroutine cududv(kidia,kfdia,klon,klev,ktopm2,ktype,kcbot,kctop,ldcum,  &
                    ptsphy,paph,puen,pven,pmfu,pmfd,puu,pud,pvu,pvd,ptenu, &
                    ptenv)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    integer(ik4) , intent(in) :: ktopm2
    integer(ik4) , dimension(klon) , intent(in) :: ktype
    integer(ik4) , dimension(klon) , intent(in) :: kcbot
    integer(ik4) , dimension(klon) , intent(in) :: kctop
    logical , dimension(klon) , intent(in) :: ldcum
    real(rk8) , intent(in) :: ptsphy
    real(rk8) , dimension(klon,klev+1) , intent(in) :: paph
    real(rk8) , dimension(klon,klev) , intent(in) :: puen
    real(rk8) , dimension(klon,klev) , intent(in) :: pven
    real(rk8) , dimension(klon,klev) , intent(in) :: pmfu
    real(rk8) , dimension(klon,klev) , intent(in) :: pmfd
    real(rk8) , dimension(klon,klev) , intent(in) :: puu
    real(rk8) , dimension(klon,klev) , intent(in) :: pud
    real(rk8) , dimension(klon,klev) , intent(in) :: pvu
    real(rk8) , dimension(klon,klev) , intent(in) :: pvd
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptenu
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptenv
    real(rk8) , dimension(klon,klev) :: zuen , zven , zmfuu , &
                 zmfdu , zmfuv , zmfdv
    integer(ik4) :: ik , ikb , jk , jl
    real(rk8) :: zzp , zimp , ztsphy
    real(rk8) , dimension(:,:) , allocatable :: zdudt , zdvdt , zdp
    real(rk8) , dimension(:,:) , allocatable :: zb , zr1 , zr2
    logical , dimension(:,:) , allocatable :: llcumbas
    zimp = d_one - rmfsoluv
    ztsphy = d_one/ptsphy
    allocate (zdudt(klon,klev))
    allocate (zdvdt(klon,klev))
    allocate (zdp(klon,klev))
    do jk = 1 , klev
      do jl = kidia , kfdia
        if ( ldcum(jl) ) then
          zuen(jl,jk) = puen(jl,jk)
          zven(jl,jk) = pven(jl,jk)
          zdp(jl,jk) = egrav/(paph(jl,jk+1)-paph(jl,jk))
        end if
      end do
    end do
    !----------------------------------------------------------------------
    !*    1.0          CALCULATE FLUXES AND UPDATE U AND V TENDENCIES
    ! ----------------------------------------------
    do jk = ktopm2 , klev
      ik = jk - 1
      do jl = kidia , kfdia
        if ( ldcum(jl) ) then
          zmfuu(jl,jk) = pmfu(jl,jk)*(puu(jl,jk)-zimp*zuen(jl,ik))
          zmfuv(jl,jk) = pmfu(jl,jk)*(pvu(jl,jk)-zimp*zven(jl,ik))
          zmfdu(jl,jk) = pmfd(jl,jk)*(pud(jl,jk)-zimp*zuen(jl,ik))
          zmfdv(jl,jk) = pmfd(jl,jk)*(pvd(jl,jk)-zimp*zven(jl,ik))
        end if
      end do
    end do
    ! linear fluxes below cloud
    if ( rmfsoluv == d_zero ) then
      do jk = ktopm2 , klev
      !DIR$ IVDEP
      !OCL NOVREC
        do jl = kidia , kfdia
          if ( ldcum(jl) .and. jk > kcbot(jl) ) then
            ikb = kcbot(jl)
            zzp = ((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
            if ( ktype(jl) == 3 ) zzp = zzp*zzp
            zmfuu(jl,jk) = zmfuu(jl,ikb)*zzp
            zmfuv(jl,jk) = zmfuv(jl,ikb)*zzp
            zmfdu(jl,jk) = zmfdu(jl,ikb)*zzp
            zmfdv(jl,jk) = zmfdv(jl,ikb)*zzp
          end if
        end do
      end do
    end if
    !*    1.2          COMPUTE TENDENCIES
    ! ------------------
    do jk = ktopm2 , klev
      if ( jk < klev ) then
        ik = jk + 1
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
    if ( rmfsoluv == d_zero ) then
      !*  1.3          UPDATE TENDENCIES
      !   -----------------
      do jk = ktopm2 , klev
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            ptenu(jl,jk) = ptenu(jl,jk) + zdudt(jl,jk)
            ptenv(jl,jk) = ptenv(jl,jk) + zdvdt(jl,jk)
          end if
        end do
      end do
    else
      !----------------------------------------------------------------------
      !*  1.6          IMPLICIT SOLUTION
      !   -----------------
      ! Fill bi-diagonal Matrix vectors A=k-1, B=k;
      ! reuse ZMFUU=A and ZB=B;
      ! ZDUDT and ZDVDT correspond to the RHS ("constants") of the equation
      ! The solution is in ZR1 and ZR2
      allocate (zb(klon,klev))
      allocate (zr1(klon,klev))
      allocate (zr2(klon,klev))
      allocate (llcumbas(klon,klev))
      llcumbas(:,:) = .false.
      zb(:,:) = d_one
      zmfuu(:,:) = d_zero
      ! Fill vectors A, B and RHS
      do jk = ktopm2 , klev
        ik = jk + 1
        do jl = kidia , kfdia
          llcumbas(jl,jk) = ldcum(jl) .and. jk >= kctop(jl) - 1
          if ( llcumbas(jl,jk) ) then
            zzp = rmfsoluv*zdp(jl,jk)*ptsphy
            zmfuu(jl,jk) = -zzp*(pmfu(jl,jk)+pmfd(jl,jk))
            zdudt(jl,jk) = zdudt(jl,jk)*ptsphy + zuen(jl,jk)
            zdvdt(jl,jk) = zdvdt(jl,jk)*ptsphy + zven(jl,jk)
            if ( jk < klev ) then
              zb(jl,jk) = d_one + zzp*(pmfu(jl,ik)+pmfd(jl,ik))
            else
              zb(jl,jk) = d_one
            end if
          end if
        end do
      end do
      call cubidiag(kidia,kfdia,klon,klev,kctop,llcumbas,zmfuu,zb,zdudt,zr1)
      call cubidiag(kidia,kfdia,klon,klev,kctop,llcumbas,zmfuu,zb,zdvdt,zr2)
      do jk = ktopm2 , klev
        do jl = kidia , kfdia
          if ( llcumbas(jl,jk) ) then
            ptenu(jl,jk) = ptenu(jl,jk) + (zr1(jl,jk)-zuen(jl,jk))*ztsphy
            ptenv(jl,jk) = ptenv(jl,jk) + (zr2(jl,jk)-zven(jl,jk))*ztsphy
          end if
        end do
      end do
      deallocate (llcumbas)
      deallocate (zr2)
      deallocate (zr1)
      deallocate (zb)
    end if
    !----------------------------------------------------------------------
    deallocate (zdp)
    deallocate (zdvdt)
    deallocate (zdudt)
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
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , dimension(klon) , intent(in) :: kctop
    logical , dimension(klon,klev) , intent(in) :: ld_lcumask
    real(rk8) , dimension(klon,klev) , intent(in) :: pa
    real(rk8) , dimension(klon,klev) , intent(in) :: pb
    real(rk8) , dimension(klon,klev) , intent(in) :: pr
    real(rk8) , dimension(klon,klev) , intent(out) :: pu
    integer(ik4) :: jk , jl
    real(rk8) :: zbet
    !----------------------------------------------------------------------
    pu(:,:) = d_zero
    ! Forward Substitution
    do jk = 2 , klev
      do jl = kidia , kfdia
        if ( ld_lcumask(jl,jk) ) then
          if ( jk == kctop(jl)-1 ) then
            zbet = d_one/(pb(jl,jk)+1.D-35)
            pu(jl,jk) = pr(jl,jk)*zbet
          else if ( jk>kctop(jl)-1 ) then
            zbet = d_one/(pb(jl,jk)+1.D-35)
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
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    integer(ik4) , intent(in) :: ktdia
    real(rk8) , dimension(klon,klev) , intent(in) :: paprsf
    real(rk8) , dimension(klon,klev) , intent(in) :: pt
    real(rk8) , dimension(klon,klev) , intent(out) :: pqsat
    integer(ik4) , intent(in) :: kflag
    integer(ik4) :: jk , jl , jlen
    real(rk8) :: zcor , zew , zqmax , zqs
    real(rk8) , dimension(kidia:kfdia) :: z_exparg1
    real(rk8) , dimension(kidia:kfdia) :: z_exparg2
    real(rk8) , dimension(kidia:kfdia) :: z_expout1
    real(rk8) , dimension(kidia:kfdia) :: z_expout2
    !----------------------------------------------------------------------
    !*    1.           DEFINE CONSTANTS
    ! ----------------
    zqmax = d_half
    ! *
    !----------------------------------------------------------------------
    ! *    2.           CALCULATE SATURATION SPECIFIC HUMIDITY
    ! --------------------------------------
    if ( n_vmass <= 0 ) then ! Not using Vector MASS
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
    else ! Using Vector MASS
      jlen = kfdia - kidia + 1
      do jk = ktdia , klev
        do jl = kidia , kfdia
          z_exparg1(jl) = foeles_v(pt(jl,jk))
          z_exparg2(jl) = foeies_v(pt(jl,jk))
        end do
        call vexp(z_expout1,z_exparg1,jlen)
        call vexp(z_expout2,z_exparg2,jlen)
        do jl = kidia , kfdia
          if ( kflag == 1 ) then
            zew = foeewmcu_v(pt(jl,jk),z_expout1(jl),z_expout2(jl))
          else
            zew = foeewm_v(pt(jl,jk),z_expout1(jl),z_expout2(jl))
          end if
          zqs = min(zqmax*paprsf(jl,jk),zew)
          pqsat(jl,jk) = zqs/(paprsf(jl,jk)-retv*zqs)
        end do
      end do
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
  subroutine cubasmcn(kidia,kfdia,klon,klev,kk,pten,pqen,pqsen,pvervel, &
                      pgeo,pgeoh,ldcum,ktype,klab,kcbot,pmfu,pmfub,     &
                      plrain,ptu,pqu,plu,pmfus,pmfuq,pmful,pdmfup)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    integer(ik4) , intent(in) :: kk
    real(rk8) , dimension(klon,klev) , intent(in) :: pten
    real(rk8) , dimension(klon,klev) , intent(in) :: pqen
    real(rk8) , dimension(klon,klev) , intent(in) :: pqsen
    real(rk8) , dimension(klon,klev) , intent(in) :: pvervel
    real(rk8) , dimension(klon,klev) , intent(in) :: pgeo
    real(rk8) , dimension(klon,klev+1) , intent(in) :: pgeoh
    logical , dimension(klon) , intent(in) :: ldcum
    integer(ik4) , dimension(klon) , intent(out) :: ktype
    integer(ik4) , dimension(klon,klev) , intent(inout) :: klab
    integer(ik4) , dimension(klon) , intent(out) :: kcbot
    real(rk8) , dimension(klon,klev) , intent(out) :: pmfu
    real(rk8) , dimension(klon) , intent(out) :: pmfub
    real(rk8) , dimension(klon,klev) , intent(out) :: plrain
    real(rk8) , dimension(klon,klev) , intent(out) :: ptu
    real(rk8) , dimension(klon,klev) , intent(out) :: pqu
    real(rk8) , dimension(klon,klev) , intent(out) :: plu
    real(rk8) , dimension(klon,klev) , intent(out) :: pmfus
    real(rk8) , dimension(klon,klev) , intent(out) :: pmfuq
    real(rk8) , dimension(klon,klev) , intent(out) :: pmful
    real(rk8) , dimension(klon,klev) , intent(out) :: pdmfup
    integer(ik4) :: jl
    real(rk8) :: zzzmb
    !----------------------------------------------------------------------
    !*    1.           CALCULATE ENTRAINMENT AND DETRAINMENT RATES
    ! -------------------------------------------
!DIR$ IVDEP
!OCL NOVREC
    do jl = kidia , kfdia
      if ( .not.ldcum(jl) .and. klab(jl,kk+1) == 0 ) then
        if ( lmfmid .and. pgeo(jl,kk) >  5000.0D0 .and. &
                          pgeo(jl,kk) < 10000.0D0 .and. &
                          pqen(jl,kk) > 0.80D0*pqsen(jl,kk) ) then
          ptu(jl,kk+1) = (rcpd*pten(jl,kk)+pgeo(jl,kk)-pgeoh(jl,kk+1))/rcpd
          pqu(jl,kk+1) = pqen(jl,kk)
          plu(jl,kk+1) = d_zero
          zzzmb = max(rmfcmin,-pvervel(jl,kk)/egrav)
          zzzmb = min(zzzmb,rmfcmax)
          pmfub(jl) = zzzmb
          pmfu(jl,kk+1) = pmfub(jl)
          pmfus(jl,kk+1) = pmfub(jl)*(rcpd*ptu(jl,kk+1)+pgeoh(jl,kk+1))
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
  subroutine cudlfsn(kidia,kfdia,klon,klev,kcbot,kctop,ldcum,ptenh,pqenh, &
                     pten,pqsen,pgeo,pgeoh,paph,ptu,pqu,pmfub,prfl,ptd,   &
                     pqd,pmfd,pmfds,pmfdq,pdmfdp,kdtop,lddraf)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    real(rk8) , dimension(klon,klev) , intent(in) :: ptenh
    real(rk8) , dimension(klon,klev) , intent(in) :: pqenh
    real(rk8) , dimension(klon,klev) , intent(in) :: pten
    real(rk8) , dimension(klon,klev) , intent(in) :: pqsen
    real(rk8) , dimension(klon,klev) , intent(in) :: pgeo
    real(rk8) , dimension(klon,klev+1) , intent(in) :: pgeoh
    real(rk8) , dimension(klon,klev+1) , intent(in) :: paph
    real(rk8) , dimension(klon,klev) , intent(in) :: ptu
    real(rk8) , dimension(klon,klev) , intent(in) :: pqu
    real(rk8) , dimension(klon) , intent(in) :: pmfub
    real(rk8) , dimension(klon) , intent(inout) :: prfl
    real(rk8) , dimension(klon,klev) , intent(out) :: ptd
    real(rk8) , dimension(klon,klev) , intent(out) :: pqd
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfd
    real(rk8) , dimension(klon,klev) , intent(out) :: pmfds
    real(rk8) , dimension(klon,klev) , intent(out) :: pmfdq
    real(rk8) , dimension(klon,klev) , intent(out) :: pdmfdp
    integer(ik4) , dimension(klon) , intent(out) :: kdtop
    integer(ik4) , dimension(klon) , intent(in) :: kctop
    integer(ik4) , dimension(klon) , intent(in) :: kcbot
    logical , dimension(klon) , intent(out) :: lddraf
    logical , dimension(klon) , intent(in) :: ldcum
    integer(ik4) , dimension(klon) :: ikhsmin
    real(rk8) , dimension(klon,klev) :: ztenwb , zqenwb
    real(rk8) , dimension(klon) :: zcond , zph , zhsmin
    logical , dimension(klon) :: llo2
    integer(ik4) :: icall , ik , ike , is , jk , jl
    real(rk8) :: zbuo , zhsk , zmftop , zqtest , zttest
    !----------------------------------------------------------------------
    ! 1.           SET DEFAULT VALUES FOR DOWNDRAFTS
    ! ---------------------------------
    do jl = kidia , kfdia
      lddraf(jl) = .false.
      kdtop(jl) = klev + 1
      ikhsmin(jl) = klev + 1
      zhsmin(jl) = 1.D8
    end do
    if ( lmfdd ) then
      !----------------------------------------------------------------------
      !   2.           DETERMINE LEVEL OF FREE SINKING:
      !   DOWNDRAFTS SHALL START AT MODEL LEVEL OF MINIMUM
      !   OF SATURATION MOIST STATIC ENERGY OR BELOW
      !   RESPECTIVELY
      !   FOR EVERY POINT AND PROCEED AS FOLLOWS:
      !   (1) DETERMINE LEVEL OF MINIMUM OF HS
      !   (2) DETERMINE WET BULB ENVIRONMENTAL T AND Q
      !   (3) DO MIXING WITH CUMULUS CLOUD AIR
      !   (4) CHECK FOR NEGATIVE BUOYANCY
      !   (5) IF BUOYANCY>0 REPEAT (2) TO (4) FOR NEXT
      !       LEVEL BELOW
      !   THE ASSUMPTION IS THAT AIR OF DOWNDRAFTS IS MIXTURE
      !   OF 50% CLOUD AIR + 50% ENVIRONMENTAL AIR AT WET BULB
      !   TEMPERATURE (I.E. WHICH BECAME SATURATED DUE TO
      !   EVAPORATION OF RAIN AND CLOUD WATER)
      !   ----------------------------------------------------
      do jk = 3 , klev - 2
        do jl = kidia , kfdia
          zhsk = rcpd*pten(jl,jk) + pgeo(jl,jk) + &
            foelhmcu(pten(jl,jk))*pqsen(jl,jk)
          if ( zhsk < zhsmin(jl) ) then
            zhsmin(jl) = zhsk
            ikhsmin(jl) = jk
          end if
        end do
      end do
      ike = klev - 3
      do jk = 3 , ike
        !     2.1          CALCULATE WET-BULB TEMPERATURE AND MOISTURE
        !     FOR ENVIRONMENTAL AIR IN *CUADJTQ*
        !     -------------------------------------------
        is = 0
        do jl = kidia , kfdia
          ztenwb(jl,jk) = ptenh(jl,jk)
          zqenwb(jl,jk) = pqenh(jl,jk)
          zph(jl) = paph(jl,jk)
          llo2(jl) = ldcum(jl) .and.                      &
            prfl(jl) > d_zero .and. .not.lddraf(jl) .and. &
            (jk < kcbot(jl) .and. jk > kctop(jl)) .and.   &
            jk >= ikhsmin(jl)
          if ( llo2(jl) ) is = is + 1
        end do
        if ( is == 0 ) cycle
        ik = jk
        icall = 2
        call cuadjtq(kidia,kfdia,klon,klev,ik,zph,ztenwb,zqenwb,llo2,icall)
        !     2.2          DO MIXING OF CUMULUS AND ENVIRONMENTAL AIR
        !     AND CHECK FOR NEGATIVE BUOYANCY.
        !     THEN SET VALUES FOR DOWNDRAFT AT LFS.
        !     ----------------------------------------
!DIR$ IVDEP
!OCL  NOVREC
        do jl = kidia , kfdia
          if ( llo2(jl) ) then
            zttest = d_half*(ptu(jl,jk)+ztenwb(jl,jk))
            zqtest = d_half*(pqu(jl,jk)+zqenwb(jl,jk))
            zbuo = zttest*(d_one+retv*zqtest) - &
                   ptenh(jl,jk)*(d_one+retv*pqenh(jl,jk))
            zcond(jl) = pqenh(jl,jk) - zqenwb(jl,jk)
            zmftop = -rmfdeps*pmfub(jl)
            if ( zbuo < d_zero .and. prfl(jl) > 10.0D0*zmftop*zcond(jl) ) then
              kdtop(jl) = jk
              lddraf(jl) = .true.
              ptd(jl,jk) = zttest
              pqd(jl,jk) = zqtest
              pmfd(jl,jk) = zmftop
              pmfds(jl,jk) = pmfd(jl,jk)*(rcpd*ptd(jl,jk)+pgeoh(jl,jk))
              pmfdq(jl,jk) = pmfd(jl,jk)*pqd(jl,jk)
              pdmfdp(jl,jk-1) = -d_half*pmfd(jl,jk)*zcond(jl)
              prfl(jl) = prfl(jl) + pdmfdp(jl,jk-1)
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
  subroutine cuddrafn(kidia,kfdia,klon,klev,lddraf,ptenh,pqenh, &
                      pgeo,pgeoh,paph,prfl,ptd,pqd,pmfu,pmfd,   &
                      pmfds,pmfdq,pdmfdp,pdmfde,pmfdde_rate,pkined)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    logical , dimension(klon) , intent(in) :: lddraf
    real(rk8) , dimension(klon,klev) , intent(in) :: ptenh
    real(rk8) , dimension(klon,klev) , intent(in) :: pqenh
    real(rk8) , dimension(klon,klev) , intent(in) :: pgeo
    real(rk8) , dimension(klon,klev+1) , intent(in) :: pgeoh
    real(rk8) , dimension(klon,klev+1) , intent(in) :: paph
    real(rk8) , dimension(klon) , intent(inout) :: prfl
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptd
    real(rk8) , dimension(klon,klev) , intent(inout) :: pqd
    real(rk8) , dimension(klon,klev) , intent(in) :: pmfu
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfd
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfds
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfdq
    real(rk8) , dimension(klon,klev) , intent(out) :: pdmfdp
    real(rk8) , dimension(klon,klev) , intent(out) :: pdmfde
    real(rk8) , dimension(klon,klev) , intent(out) :: pmfdde_rate
    real(rk8) , dimension(klon,klev) , intent(out) :: pkined
    real(rk8) , dimension(klon) :: zdmfen , zdmfde , zcond , zoentr , zbuoy
    real(rk8) , dimension(klon) :: zph
    logical , dimension(klon) :: llo2
    integer(ik4) :: icall , ik , is , itopde , jk , jl
    real(rk8) :: zbuo , zbuoyz , zbuoyv , zdmfdp , zdz , zentr , zmfdqk ,  &
                 zmfdsk , zqdde , zqeen , zrain , zsdde , zseen , zzentr , &
                 zfacbuo , z_cwdrag , zdkbuo , zdken
    itopde = njkt3
    zfacbuo = d_half/(d_one+d_half)
    z_cwdrag = (3.0D0/8.0D0)*0.506D0/0.200D0
    !----------------------------------------------------------------------
    ! 1.           CALCULATE MOIST DESCENT FOR CUMULUS DOWNDRAFT BY
    ! (A) CALCULATING ENTRAINMENT/DETRAINMENT RATES,
    ! INCLUDING ORGANIZED ENTRAINMENT DEPENDENT ON
    ! NEGATIVE BUOYANCY AND ASSUMING
    ! LINEAR DECREASE OF MASSFLUX IN PBL
    ! (B) DOING MOIST DESCENT - EVAPORATIVE COOLING
    ! AND MOISTENING IS CALCULATED IN *CUADJTQ*
    ! (C) CHECKING FOR NEGATIVE BUOYANCY AND
    ! SPECIFYING FINAL T,Q,U,V AND DOWNWARD FLUXES
    ! -------------------------------------------------
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
        if ( llo2(jl) ) is = is + 1
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
            zdmfen(jl) = zdmfen(jl) + zzentr
            zdmfen(jl) = max(zdmfen(jl),0.3D0*pmfd(jl,jk-1))
            zdmfen(jl) = max(zdmfen(jl), &
              -0.75D0*pmfu(jl,jk)-(pmfd(jl,jk-1)-zdmfde(jl)))
            zdmfen(jl) = min(zdmfen(jl),d_zero)
          end if
          pdmfde(jl,jk) = zdmfen(jl) - zdmfde(jl)
        end do
      end if
      do jl = kidia , kfdia
        if ( llo2(jl) ) then
          pmfd(jl,jk) = pmfd(jl,jk-1) + zdmfen(jl) - zdmfde(jl)
          zseen = (rcpd*ptenh(jl,jk-1)+pgeoh(jl,jk-1))*zdmfen(jl)
          zqeen = pqenh(jl,jk-1)*zdmfen(jl)
          zsdde = (rcpd*ptd(jl,jk-1)+pgeoh(jl,jk-1))*zdmfde(jl)
          zqdde = pqd(jl,jk-1)*zdmfde(jl)
          zmfdsk = pmfds(jl,jk-1) + zseen - zsdde
          zmfdqk = pmfdq(jl,jk-1) + zqeen - zqdde
          pqd(jl,jk) = zmfdqk*(d_one/min(-rmfcmin,pmfd(jl,jk)))
          ptd(jl,jk) = (zmfdsk*(d_one / &
            min(-rmfcmin,pmfd(jl,jk)))-pgeoh(jl,jk))/rcpd
          ptd(jl,jk) = min(400.0D0,ptd(jl,jk))
          ptd(jl,jk) = max(100.0D0,ptd(jl,jk))
          zcond(jl) = pqd(jl,jk)
        end if
      end do
      ik = jk
      icall = 2
      call cuadjtq(kidia,kfdia,klon,klev,ik,zph,ptd,pqd,llo2,icall)
      do jl = kidia , kfdia
        if ( llo2(jl) ) then
          zcond(jl) = zcond(jl) - pqd(jl,jk)
          zbuo = ptd(jl,jk)*(d_one+retv*pqd(jl,jk)) - &
            ptenh(jl,jk)*(d_one+retv*pqenh(jl,jk))
          if ( prfl(jl) > d_zero .and. pmfu(jl,jk) > d_zero ) then
            zrain = prfl(jl)/pmfu(jl,jk)
            zbuo = zbuo - ptd(jl,jk)*zrain
          end if
          if ( zbuo >= d_zero .or. prfl(jl)<=(pmfd(jl,jk)*zcond(jl)) ) then
            pmfd(jl,jk) = d_zero
            zbuo = d_zero
          end if
          pmfds(jl,jk) = (rcpd*ptd(jl,jk)+pgeoh(jl,jk))*pmfd(jl,jk)
          pmfdq(jl,jk) = pqd(jl,jk)*pmfd(jl,jk)
          zdmfdp = -pmfd(jl,jk)*zcond(jl)
          pdmfdp(jl,jk-1) = zdmfdp
          prfl(jl) = prfl(jl) + zdmfdp
          !       COMPUTE ORGANIZED ENTRAINMENT FOR USE AT NEXT LEVEL
          zbuoyz = zbuo/ptenh(jl,jk)
          zbuoyv = zbuoyz
          zbuoyz = min(zbuoyz,d_zero)
          zdz = -(pgeo(jl,jk-1)-pgeo(jl,jk))
          zbuoy(jl) = zbuoy(jl) + zbuoyz*zdz
          zoentr(jl) = egrav*zbuoyz*d_half/(d_one+zbuoy(jl))
          !       STORE DOWNDRAUGHT DETRAINMENT RATES
          pmfdde_rate(jl,jk) = -zdmfde(jl)
          !       COMPUTE KINETIC ENERGY
          zdkbuo = zdz*zbuoyv*zfacbuo
          if ( zdmfen(jl) < d_zero ) then
            zdken = min(d_one,(d_one+z_cwdrag)*zdmfen(jl) / &
                    min(-rmfcmin,pmfd(jl,jk-1)))
          else
            zdken = min(d_one,(d_one+z_cwdrag)*zdmfde(jl) / &
                    min(-rmfcmin,pmfd(jl,jk-1)))
          end if
          pkined(jl,jk) = max(d_zero, &
            (pkined(jl,jk-1)*(d_one-zdken)+zdkbuo)/(d_one+zdken))
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
  subroutine custrat(kidia,kfdia,klon,klev,ldcum,ptsphy,pap,paph,pgeo,   &
                     pten,pqen,pqsat,penth,ptent,ptenq)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    logical , dimension(klon) , intent(in) :: ldcum(klon)
    real(rk8) , intent(in) :: ptsphy
    real(rk8) , dimension(klon,klev) , intent(in) :: pap
    real(rk8) , dimension(klon,klev+1) , intent(in) :: paph
    real(rk8) , dimension(klon,klev) , intent(in) :: pgeo
    real(rk8) , dimension(klon,klev) , intent(in) :: pten
    real(rk8) , dimension(klon,klev) , intent(in) :: pqen
    real(rk8) , dimension(klon,klev) , intent(in) :: pqsat
    real(rk8) , dimension(klon,klev) , intent(out) :: penth
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptent
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptenq
    real(rk8) , dimension(klon,klev) :: ztc , zqc , zcf, zcptgz , &
      ztdif , zqdif , zebs
    real(rk8) , dimension(klon,klev+1) :: zap
    real(rk8) , dimension(klon) :: zcpts , zqs , ztcoe , zqold
    real(rk8) , dimension(klon) :: zpp
    integer(ik4) , dimension(klon,klev) :: ilab
    logical , dimension(klon) :: llflag , llo2 , llbl
    integer(ik4) :: icall , ik , ilevh , jk , jl
    real(rk8) :: zbuo , zcons1 , zcons2 , zcons3 , zdisc , zdqdt , zdtdt , &
                 zfac , zkdiff1 , zkdiff2 , zqdp , ztmst , ztpfac1 , ztpfac2
    !-----------------------------------------------------------------------
    ! 2.           PHYSICAL CONSTANTS AND PARAMETERS.
    ! ---------------------------------
    ztpfac1 = rvdifts
    ztpfac2 = d_one/ztpfac1
    zkdiff1 = 10.0D0
    zkdiff2 = 2.5D0
    ztmst = ptsphy
    zcons1 = ztpfac1*ztmst*egrav**2/(d_half*rgas)
    zcons2 = d_one/ztmst
    zcons3 = ztmst*rcpd
    ilevh = klev/2
    !----------------------------------------------------------------------
    !*    3.           PRELIMINARY COMPUTATIONS.
    ! ------------------------
    do jk = 1 , klev
      do jl = kidia , kfdia
        zcptgz(jl,jk) = pgeo(jl,jk) + pten(jl,jk)*rcpd
        zcf(jl,jk) = d_zero
        ilab(jl,jk) = 0
      end do
    end do
    !-----------------------------------------------------------------
    ! 4.           DETERMINE EXCHANGE COEFFICIENTS THEREFORE
    ! (A) LIFT SURFACE AIR, CHECK FOR BUOYANCY AND SET FLAG
    ! (B) THEN DEFINE DIFFUSION COEFFICIENTS,I.E.
    ! K=C1 FOR CLOUD LAYER
    ! K=C1*F(RH) FOR CLOUD TOP (TOP ENTRAINMENT)
    ! ----------------------------------------------------
    do jl = kidia , kfdia
      ztc(jl,klev) = pten(jl,klev) + 0.25D0
      zqc(jl,klev) = pqen(jl,klev)
      if ( .not.ldcum(jl) ) then
        ilab(jl,klev) = 1
      else
        ilab(jl,klev) = 0
      end if
      llo2(jl) = .false.
      llbl(jl) = .true.
    end do
    do jk = klev - 1 , ilevh , -1
      do jl = kidia , kfdia
        if ( pap(jl,jk) < 0.9D0*paph(jl,klev+1) ) llbl(jl) = .false.
      end do
      do jl = kidia , kfdia
        if ( llbl(jl) ) then
          ztc(jl,jk) = (ztc(jl,jk+1)*rcpd+pgeo(jl,jk+1)-pgeo(jl,jk))/rcpd
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
      call cuadjtq(kidia,kfdia,klon,klev,ik,zpp,ztc,zqc,llflag,icall)
      do jl = kidia , kfdia
        if ( llbl(jl) ) then
          if ( zqc(jl,jk) /= zqold(jl) ) ilab(jl,jk) = 2
        end if
      end do
!DIR$ IVDEP
!OCL NOVREC
      do jl = kidia , kfdia
        if ( llbl(jl) ) then
          zbuo = ztc(jl,jk)*(d_one+retv*zqc(jl,jk)) - &
            pten(jl,jk)*(d_one+retv*pqen(jl,jk))
          if ( zbuo < d_zero ) ilab(jl,jk) = 0
          if ( zbuo > d_zero .and. &
            ilab(jl,jk) == 0 .and. ilab(jl,jk+1) == 1 ) ilab(jl,jk) = 1
          if ( ilab(jl,jk) == 2 ) llo2(jl) = .true.
        end if
      end do
    end do
    do jl = kidia , kfdia
      llbl(jl) = .true.
    end do
    do jk = klev - 1 , ilevh , -1
      do jl = kidia , kfdia
        if ( pap(jl,jk) < 0.9D0*paph(jl,klev+1) ) llbl(jl) = .false.
      end do
      do jl = kidia , kfdia
        if ( llbl(jl) ) then
          if ( ilab(jl,jk) == 2 ) then
            zcf(jl,jk) = zkdiff1
            if ( ilab(jl,klev-2) == 0 ) zcf(jl,jk) = zkdiff2
          else
            zcf(jl,jk) = d_zero
          end if
          if ( zcf(jl,jk+1) > d_zero .and. ilab(jl,jk) == 0 ) then
            zcf(jl,jk) = zcf(jl,jk+1) * &
              5.0D0*max(pqen(jl,jk+1)/pqsat(jl,jk+1)-0.8D0,d_zero) *&
                    max(pqen(jl,jk+1)/pqsat(jl,jk+1)- &
                         pqen(jl,jk)/pqsat(jl,jk),d_zero)
            llbl(jl) = .false.
          end if
        end if
      end do
    end do
    !*    4.7          EXCHANGE COEFFICIENTS.
    do jk = ilevh , klev - 1
      do jl = kidia , kfdia
        zcf(jl,jk) = zcf(jl,jk)*zcons1*paph(jl,jk+1) / &
                     ((pgeo(jl,jk)-pgeo(jl,jk+1))*(pten(jl,jk)+pten(jl,jk+1)))
      end do
    end do
    !*    4.8          DUMMY SURFACE VALUES OF T AND Q AT SURFACE
    do jl = kidia , kfdia
      zcpts(jl) = ztpfac2*zcptgz(jl,klev)
      zqs(jl) = ztpfac2*pqen(jl,klev)
    end do
    !----------------------------------------------------------------------
    ! 5.           SOLUTION OF THE VERTICAL DIFFUSION EQUATION.
    ! --------------------------------------------
    !*    5.1          SETTING OF RIGHT HAND SIDES.
    do jk = ilevh , klev
      do jl = kidia , kfdia
        ztdif(jl,jk) = ztpfac2*zcptgz(jl,jk)
        zqdif(jl,jk) = ztpfac2*pqen(jl,jk)
      end do
    end do
    !*    5.2          TOP LAYER ELIMINATION.
    do jl = kidia , kfdia
      ztcoe(jl) = zcf(jl,ilevh)
      zqdp = d_one/(paph(jl,ilevh+1)-paph(jl,ilevh))
      zdisc = d_one/(d_one+zcf(jl,ilevh)*zqdp)
      zebs(jl,ilevh) = zdisc*(zcf(jl,ilevh)*zqdp)
      zqdif(jl,ilevh) = zdisc*zqdif(jl,ilevh)
      ztdif(jl,ilevh) = zdisc*ztdif(jl,ilevh)
    end do
    !*    5.3          ELIMINATION FOR LAYERS BELOW
    do jk = ilevh + 1 , klev
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
      zqdif(jl,klev) = zqdif(jl,klev) + (zebs(jl,klev)*zqs(jl))
      ztdif(jl,klev) = ztdif(jl,klev) + (zebs(jl,klev)*zcpts(jl))
    end do
    !*    5.5          BACK-SUBSTITUTION.
    do jk = klev - 1 , ilevh , -1
      do jl = kidia , kfdia
        zqdif(jl,jk) = zqdif(jl,jk) + (zebs(jl,jk)*zqdif(jl,jk+1))
        ztdif(jl,jk) = ztdif(jl,jk) + (zebs(jl,jk)*ztdif(jl,jk+1))
      end do
    end do
    !---------------------------------------------------------------------
    !*    6.           INCREMENTATION OF T AND Q TENDENCIES.
    ! -------------------------------------
    do jk = ilevh , klev
      do jl = kidia , kfdia
        zdqdt = (zqdif(jl,jk)-ztpfac2*pqen(jl,jk))*zcons2
        ptenq(jl,jk) = ptenq(jl,jk) + zdqdt
        zdtdt = (ztdif(jl,jk)-ztpfac2*zcptgz(jl,jk))/zcons3
        ptent(jl,jk) = ptent(jl,jk) + zdtdt
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
  subroutine cudtdqn(kidia,kfdia,klon,klev,ktopm2,kctop,kdtop,ldcum, &
                     lddraf,ptsphy,paph,pgeoh,pgeo,pten,ptenh,pqen,  &
                     pqenh,pqsen,plglac,plude,pmfu,pmfd,pmfus,pmfds, &
                     pmfuq,pmfdq,pmful,pdmfup,pdpmel,ptent,ptenq,penth)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    integer(ik4) , intent(in) :: ktopm2
    integer(ik4) , dimension(klon) , intent(in) :: kctop
    integer(ik4) , dimension(klon) , intent(in) :: kdtop
    logical , dimension(klon) , intent(inout) :: ldcum
    logical , dimension(klon) , intent(in) :: lddraf
    real(rk8) , intent(in) :: ptsphy
    real(rk8) , dimension(klon,klev+1) , intent(in) :: paph
    real(rk8) , dimension(klon,klev) , intent(in) :: pgeo
    real(rk8) , dimension(klon,klev+1) , intent(in) :: pgeoh
    real(rk8) , dimension(klon,klev) , intent(in) :: pten
    real(rk8) , dimension(klon,klev) , intent(in) :: pqen
    real(rk8) , dimension(klon,klev) , intent(in) :: ptenh
    real(rk8) , dimension(klon,klev) , intent(in) :: pqenh
    real(rk8) , dimension(klon,klev) , intent(in) :: pqsen
    real(rk8) , dimension(klon,klev) , intent(in) :: plglac
    real(rk8) , dimension(klon,klev) , intent(inout) :: plude
    real(rk8) , dimension(klon,klev) , intent(in) :: pmfu
    real(rk8) , dimension(klon,klev) , intent(in) :: pmfd
    real(rk8) , dimension(klon,klev) , intent(in) :: pmfus
    real(rk8) , dimension(klon,klev) , intent(in) :: pmfds
    real(rk8) , dimension(klon,klev) , intent(in) :: pmfuq
    real(rk8) , dimension(klon,klev) , intent(in) :: pmfdq
    real(rk8) , dimension(klon,klev) , intent(in) :: pmful
    real(rk8) , dimension(klon,klev) , intent(in) :: pdmfup
    real(rk8) , dimension(klon,klev) , intent(in) :: pdpmel
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptent
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptenq
    real(rk8) , dimension(klon,klev) , intent(out) :: penth
    logical :: lltest
    integer(ik4) :: jk , ik , jl
    real(rk8) :: ztsphy , zimp , zalv , zzp , zgq , zgs , zgh , zs , zq
    real(rk8) , dimension(klon,klev) :: zmfus , zmfuq , zmfds , zmfdq
    real(rk8) , dimension(:,:) , allocatable :: zdtdt , zdqdt , zdp
    real(rk8) , dimension(:,:) , allocatable :: zb , zr1 , zr2
    logical , dimension(:,:) , allocatable :: llcumbas
    !*    1.0          SETUP AND INITIALIZATIONS
    ! -------------------------
    zimp = d_one - rmfsoltq
    ztsphy = d_one/ptsphy
    allocate (zdtdt(klon,klev))
    allocate (zdqdt(klon,klev))
    allocate (zdp(klon,klev))
    do jk = 1 , klev
      do jl = kidia , kfdia
        penth(jl,jk) = d_zero
      end do
    end do
    ! zero detrained liquid water if diagnostic cloud scheme to be used
    ! this means that detrained liquid water will be evaporated in the
    ! cloud environment and not fed directly into a cloud liquid water
    ! variable
    lltest = .not. lepcld
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
    !-----------------------------------------------------------------------
    if ( rmfsoltq > d_zero ) then
      !*  2.0          RECOMPUTE CONVECTIVE FLUXES IF IMPLICIT
      do jk = ktopm2 , klev
        ik = jk - 1
!DIR$ IVDEP
!OCL NOVREC
        do jl = kidia , kfdia
          if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
            ! compute interpolating coefficients ZGS and ZGQ
            ! for half-level values
            zgq = (pqenh(jl,jk)-pqen(jl,ik))/pqsen(jl,jk)
            zgh = rcpd*pten(jl,jk) + pgeo(jl,jk)
            zgs = (rcpd*(ptenh(jl,jk)-pten(jl,ik)) + &
              pgeoh(jl,jk)-pgeo(jl,ik))/zgh
            !half-level environmental values for S and Q
            zs = rcpd*(zimp*pten(jl,ik)+zgs*pten(jl,jk)) + &
              pgeo(jl,ik) + zgs*pgeo(jl,jk)
            zq = zimp*pqen(jl,ik) + zgq*pqsen(jl,jk)
            zmfus(jl,jk) = pmfus(jl,jk) - pmfu(jl,jk)*zs
            zmfuq(jl,jk) = pmfuq(jl,jk) - pmfu(jl,jk)*zq
            if ( lddraf(jl) .and. jk >= kdtop(jl) ) then
              zmfds(jl,jk) = pmfds(jl,jk) - pmfd(jl,jk)*zs
              zmfdq(jl,jk) = pmfdq(jl,jk) - pmfd(jl,jk)*zq
            end if
          end if
        end do
      end do
    end if
    !*    3.0          COMPUTE TENDENCIES
    ! ------------------
    do jk = ktopm2 , klev
      if ( jk < klev ) then
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            zalv = foelhmcu(pten(jl,jk))
            zdtdt(jl,jk) = zdp(jl,jk)*rcpd * &
              (zmfus(jl,jk+1)-zmfus(jl,jk)+zmfds(jl,jk+1) - &
               zmfds(jl,jk)+wlhf*plglac(jl,jk)-wlhf*pdpmel(jl,jk) - &
               zalv*(pmful(jl,jk+1)-pmful(jl,jk)-plude(jl,jk)-pdmfup(jl,jk)))
            zdqdt(jl,jk) = zdp(jl,jk)*(zmfuq(jl,jk+1) - &
              zmfuq(jl,jk)+zmfdq(jl,jk+1)-zmfdq(jl,jk)+pmful(jl,jk+1) - &
              pmful(jl,jk)-plude(jl,jk)-pdmfup(jl,jk))
          end if
        end do
      else
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            zalv = foelhmcu(pten(jl,jk))
            zdtdt(jl,jk) = -zdp(jl,jk)*rcpd * &
              (zmfus(jl,jk)+zmfds(jl,jk)+wlhf*pdpmel(jl,jk) - &
               zalv*(pmful(jl,jk)+pdmfup(jl,jk)))
            zdqdt(jl,jk) = -zdp(jl,jk)*(zmfuq(jl,jk) + &
              zmfdq(jl,jk)+(pmful(jl,jk)+pdmfup(jl,jk)))
          end if
        end do
      end if
    end do
    if ( rmfsoltq == d_zero ) then
      !*  3.1          UPDATE TENDENCIES
      !   -----------------
      do jk = ktopm2 , klev
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            ptent(jl,jk) = ptent(jl,jk) + zdtdt(jl,jk)
            ptenq(jl,jk) = ptenq(jl,jk) + zdqdt(jl,jk)
            penth(jl,jk) = zdtdt(jl,jk)*rcpd
          end if
        end do
      end do
    else
      !----------------------------------------------------------------------
      !*  3.2          IMPLICIT SOLUTION
      !   -----------------
      ! Fill bi-diagonal Matrix vectors A=k-1, B=k, C=k+1;
      ! reuse ZMFUS=A
      ! ZDTDT and ZDQDT correspond to the RHS ("constants") of the equation
      ! The solution is in ZR1 and ZR2
      allocate (zb(klon,klev))
      allocate (zr1(klon,klev))
      allocate (zr2(klon,klev))
      allocate (llcumbas(klon,klev))
      llcumbas(:,:) = .false.
      zb(:,:) = d_one
      zmfus(:,:) = d_zero
      ! Fill vectors A, B and RHS
      do jk = ktopm2 , klev
        ik = jk + 1
        do jl = kidia , kfdia
          llcumbas(jl,jk) = ldcum(jl) .and. jk >= kctop(jl) - 1
          if ( llcumbas(jl,jk) ) then
            zzp = rmfsoltq*zdp(jl,jk)*ptsphy
            zmfus(jl,jk) = -zzp*(pmfu(jl,jk)+pmfd(jl,jk))
            zdtdt(jl,jk) = zdtdt(jl,jk)*ptsphy + pten(jl,jk)
            zdqdt(jl,jk) = zdqdt(jl,jk)*ptsphy + pqen(jl,jk)
            if ( jk < klev ) then
              zb(jl,jk) = d_one + zzp*(pmfu(jl,ik)+pmfd(jl,ik))
            else
              zb(jl,jk) = d_one
            end if
          end if
        end do
      end do
      call cubidiag(kidia,kfdia,klon,klev,kctop,llcumbas,zmfus,zb,zdtdt,zr1)
      call cubidiag(kidia,kfdia,klon,klev,kctop,llcumbas,zmfus,zb,zdqdt,zr2)
      ! Compute tendencies
      do jk = ktopm2 , klev
        do jl = kidia , kfdia
          if ( llcumbas(jl,jk) ) then
            ptent(jl,jk) = ptent(jl,jk) + (zr1(jl,jk)-pten(jl,jk))*ztsphy
            ptenq(jl,jk) = ptenq(jl,jk) + (zr2(jl,jk)-pqen(jl,jk))*ztsphy
            penth(jl,jk) = (zr1(jl,jk)-pten(jl,jk))*ztsphy
          end if
        end do
      end do
      deallocate (llcumbas)
      deallocate (zr2)
      deallocate (zr1)
      deallocate (zb)
      !----------------------------------------------------------------------
    end if
    deallocate (zdp)
    deallocate (zdqdt)
    deallocate (zdtdt)
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
  subroutine cuflxn(kidia,kfdia,klon,klev,ptsphy,pten,pqen,pqsen,ptenh,  &
                    pqenh,paph,pap,pgeoh,ldland,ldcum,kcbot,kctop,kdtop, &
                    ktopm2,ktype,lddraf,pmfu,pmfd,pmfus,pmfds,pmfuq,     &
                    pmfdq,pmful,plude,pdmfup,pdmfdp,pdpmel,plglac,pmflxr,&
                    pmflxs,prain,pmfdde_rate)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    real(rk8) , intent(in) :: ptsphy
    real(rk8) , dimension(klon,klev) , intent(in) :: pten
    real(rk8) , dimension(klon,klev) , intent(in) :: pqen
    real(rk8) , dimension(klon,klev) , intent(inout) :: pqsen
    real(rk8) , dimension(klon,klev) , intent(in) :: ptenh
    real(rk8) , dimension(klon,klev) , intent(in) :: pqenh
    real(rk8) , dimension(klon,klev+1) , intent(in) :: paph
    real(rk8) , dimension(klon,klev) , intent(in) :: pap
    real(rk8) , dimension(klon,klev+1) , intent(in) :: pgeoh
    logical , dimension(klon) , intent(in) :: ldland
    logical , dimension(klon) , intent(in) :: ldcum
    integer(ik4) , dimension(klon) , intent(in) :: kcbot
    integer(ik4) , dimension(klon) , intent(in) :: kctop
    integer(ik4) , dimension(klon) , intent(in) :: kdtop
    integer(ik4) , intent(out) :: ktopm2
    integer(ik4) , dimension(klon) , intent(inout) :: ktype
    logical , dimension(klon) , intent(inout) :: lddraf
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfu
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfd
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfus
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfds
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfuq
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfdq
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmful
    real(rk8) , dimension(klon,klev) , intent(out) :: plude
    real(rk8) , dimension(klon,klev) , intent(inout) :: pdmfup
    real(rk8) , dimension(klon,klev) , intent(inout) :: pdmfdp
    real(rk8) , dimension(klon,klev) , intent(out) :: pdpmel
    real(rk8) , dimension(klon,klev) , intent(inout) :: plglac
    real(rk8) , dimension(klon,klev+1) , intent(out) :: pmflxr
    real(rk8) , dimension(klon,klev+1) , intent(out) :: pmflxs
    real(rk8) , dimension(klon) , intent(out) :: prain
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfdde_rate
    real(rk8) , dimension(klon) :: zrhebc
    integer(ik4) :: ik , ikb , jk , jl
    integer(ik4) , dimension(klon) :: idbas
    logical :: llddraf
    real(rk8) :: zalfaw , zcons1 , zcons1a , zcons2 , zdenom , zdrfl , &
                 zdrfl1 , zfac , zpdr , zpds , zrfl , zrfln , zrmin ,  &
                 zrnew , zsnmlt , ztmst , zzp
    ztmst = ptsphy
    zcons1a = rcpd/(wlhf*egrav*rtaumel)
    zcons2 = rmfcfl/(egrav*ztmst)
    !*    1.0          DETERMINE FINAL CONVECTIVE FLUXES
    ! ---------------------------------
    do jl = kidia , kfdia
      prain(jl) = d_zero
      if ( .not.ldcum(jl) .or. kdtop(jl) < kctop(jl) ) lddraf(jl) = .false.
      if ( .not.ldcum(jl) ) ktype(jl) = 0
      idbas(jl) = klev
      if ( ldland(jl) ) then
        zrhebc(jl) = 0.7D0
      else
        zrhebc(jl) = 0.9D0
      end if
    end do
    ! TO GET IDENTICAL RESULTS FOR DIFFERENT NPROMA FORCE KTOPM2 TO 2
    ktopm2 = 2
    do jk = ktopm2 , klev
!DIR$ IVDEP
!OCL NOVREC
      ikb = min(jk+1,klev)
      do jl = kidia , kfdia
        pmflxr(jl,jk) = d_zero
        pmflxs(jl,jk) = d_zero
        pdpmel(jl,jk) = d_zero
        if ( ldcum(jl) .and. jk >= kctop(jl) ) then
          pmfus(jl,jk) = pmfus(jl,jk) - &
            pmfu(jl,jk)*(rcpd*ptenh(jl,jk)+pgeoh(jl,jk))
          pmfuq(jl,jk) = pmfuq(jl,jk) - pmfu(jl,jk)*pqenh(jl,jk)
          plglac(jl,jk) = pmfu(jl,jk)*plglac(jl,jk)
          llddraf = lddraf(jl) .and. jk >= kdtop(jl)
          if ( llddraf ) then
            pmfds(jl,jk) = pmfds(jl,jk) - &
              pmfd(jl,jk)*(rcpd*ptenh(jl,jk)+pgeoh(jl,jk))
            pmfdq(jl,jk) = pmfdq(jl,jk) - pmfd(jl,jk)*pqenh(jl,jk)
          else
            pmfd(jl,jk) = d_zero
            pmfds(jl,jk) = d_zero
            pmfdq(jl,jk) = d_zero
            pdmfdp(jl,jk-1) = d_zero
          end if
          if ( llddraf .and. &
               pmfd(jl,jk) < d_zero .and. pmfd(jl,ikb) == d_zero) then
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
    !*    1.5          SCALE FLUXES BELOW CLOUD BASE
    ! LINEAR DCREASE
    ! -----------------------------
!DIR$ IVDEP
!OCL NOVREC
    do jl = kidia , kfdia
      if ( ldcum(jl) ) then
        ikb = kcbot(jl)
        ik = ikb + 1
        zzp = ((paph(jl,klev+1)-paph(jl,ik))/(paph(jl,klev+1)-paph(jl,ikb)))
        if ( ktype(jl) == 3 ) zzp = zzp*zzp
        pmfu(jl,ik) = pmfu(jl,ikb)*zzp
        pmfus(jl,ik) = (pmfus(jl,ikb)-foelhmcu(ptenh(jl,ikb))*pmful(jl,ikb))*zzp
        pmfuq(jl,ik) = (pmfuq(jl,ikb)+pmful(jl,ikb))*zzp
        pmful(jl,ik) = d_zero
      end if
    end do
    do jk = ktopm2 , klev
!DIR$ IVDEP
!OCL NOVREC
      do jl = kidia , kfdia
        if ( ldcum(jl) .and. jk > kcbot(jl)+1 ) then
          ikb = kcbot(jl) + 1
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
    !*    2.            CALCULATE RAIN/SNOW FALL RATES
    !*                  CALCULATE MELTING OF SNOW
    !*                  CALCULATE EVAPORATION OF PRECIP
    ! -------------------------------
    do jk = ktopm2 , klev
      do jl = kidia , kfdia
        if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
          prain(jl) = prain(jl) + pdmfup(jl,jk)
          if ( pmflxs(jl,jk) > d_zero .and. pten(jl,jk) > tzero ) then
            zcons1 = zcons1a*(d_one+d_half*(pten(jl,jk)-tzero))
            zfac = zcons1*(paph(jl,jk+1)-paph(jl,jk))
            zsnmlt = min(pmflxs(jl,jk),zfac*(pten(jl,jk)-tzero))
            pdpmel(jl,jk) = zsnmlt
            pqsen(jl,jk) = foeewmcu(pten(jl,jk)-zsnmlt/zfac)/pap(jl,jk)
          end if
          zalfaw = foealfcu(pten(jl,jk))
          ! no liquid precipitation above melting level
          if ( pten(jl,jk) < tzero .and. zalfaw > d_zero ) then
            plglac(jl,jk) = plglac(jl,jk) + zalfaw*(pdmfup(jl,jk)+pdmfdp(jl,jk))
            zalfaw = d_zero
          end if
          pmflxr(jl,jk+1) = pmflxr(jl,jk) + &
            zalfaw*(pdmfup(jl,jk)+pdmfdp(jl,jk))+pdpmel(jl,jk)
          pmflxs(jl,jk+1) = pmflxs(jl,jk) + &
            (d_one-zalfaw)*(pdmfup(jl,jk)+pdmfdp(jl,jk)) - pdpmel(jl,jk)
          if ( pmflxr(jl,jk+1)+pmflxs(jl,jk+1) < d_zero ) then
            pdmfdp(jl,jk) = -(pmflxr(jl,jk)+pmflxs(jl,jk)+pdmfup(jl,jk))
            pmflxr(jl,jk+1) = d_zero
            pmflxs(jl,jk+1) = d_zero
            pdpmel(jl,jk) = d_zero
          else if ( pmflxr(jl,jk+1) < d_zero ) then
            pmflxs(jl,jk+1) = pmflxs(jl,jk+1) + pmflxr(jl,jk+1)
            pmflxr(jl,jk+1) = d_zero
          else if ( pmflxs(jl,jk+1) < d_zero ) then
            pmflxr(jl,jk+1) = pmflxr(jl,jk+1) + pmflxs(jl,jk+1)
            pmflxs(jl,jk+1) = d_zero
          end if
        end if
      end do
    end do
    do jk = ktopm2 , klev
      do jl = kidia , kfdia
        if ( ldcum(jl) .and. jk >= kcbot(jl) ) then
          zrfl = pmflxr(jl,jk) + pmflxs(jl,jk)
          if ( zrfl > 1.D-20 ) then
            zdrfl1 = rcpecons * &
              max(d_zero,pqsen(jl,jk)-pqen(jl,jk))*rcucov * &
              (sqrt(paph(jl,jk)/paph(jl,klev+1)) / &
                    5.09D-3*zrfl/rcucov)**0.5777D0*(paph(jl,jk+1)-paph(jl,jk))
            zrnew = zrfl - zdrfl1
            zrmin = zrfl - rcucov*max(d_zero,zrhebc(jl)*pqsen(jl,jk) - &
                    pqen(jl,jk))*zcons2*(paph(jl,jk+1)-paph(jl,jk))
            zrnew = max(zrnew,zrmin)
            zrfln = max(zrnew,d_zero)
            zdrfl = min(d_zero,zrfln-zrfl)
            zalfaw = foealfcu(pten(jl,jk))
            if ( pten(jl,jk) < tzero ) zalfaw = d_zero
            zpdr = zalfaw*pdmfdp(jl,jk)
            zpds = (d_one-zalfaw)*pdmfdp(jl,jk)
            zdenom = d_one/max(1.D-20,pmflxr(jl,jk)+pmflxs(jl,jk))
            pmflxr(jl,jk+1) = pmflxr(jl,jk) + zpdr + pdpmel(jl,jk) + &
                              zdrfl*pmflxr(jl,jk)*zdenom
            pmflxs(jl,jk+1) = pmflxs(jl,jk) + zpds - pdpmel(jl,jk) + &
                              zdrfl*pmflxs(jl,jk)*zdenom
            pdmfup(jl,jk) = pdmfup(jl,jk) + zdrfl
            if ( pmflxr(jl,jk+1)+pmflxs(jl,jk+1) < d_zero ) then
              pdmfup(jl,jk) = pdmfup(jl,jk) - (pmflxr(jl,jk+1)+pmflxs(jl,jk+1))
              pmflxr(jl,jk+1) = d_zero
              pmflxs(jl,jk+1) = d_zero
              pdpmel(jl,jk) = d_zero
            else if ( pmflxr(jl,jk+1) < d_zero ) then
              pmflxs(jl,jk+1) = pmflxs(jl,jk+1) + pmflxr(jl,jk+1)
              pmflxr(jl,jk+1) = d_zero
            else if ( pmflxs(jl,jk+1) < d_zero ) then
              pmflxr(jl,jk+1) = pmflxr(jl,jk+1) + pmflxs(jl,jk+1)
              pmflxs(jl,jk+1) = d_zero
            end if
          else
            pmflxr(jl,jk+1) = d_zero
            pmflxs(jl,jk+1) = d_zero
            pdmfdp(jl,jk) = d_zero
            pdpmel(jl,jk) = d_zero
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

    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    integer(ik4) :: ktdia ! argument not used
    logical , dimension(klon) , intent(in) :: ldland
    real(rk8) , dimension(klon,klev) , intent(in) :: ptenh , pqenh , &
              pten , pqen , pqsen , pgeo , puen , pven
    real(rk8) , dimension(klon,klev+1) , intent(in) :: pgeoh , paph , &
              pqhfl , pahfs
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptu , pqu , &
              plu , puu , pvu
    real(rk8) , dimension(klon) , intent(out) :: pwubase , pcape
    integer(ik4) , dimension(klon,klev) , intent(inout) :: klab
    logical , dimension(klon) , intent(inout) :: ldcum
    logical , dimension(klon) , intent(out) :: ldsc
    integer(ik4) , dimension(klon) , intent(inout) :: kcbot
    integer(ik4) , dimension(klon) , intent(out) :: kbotsc , kctop , kdpl

    integer(ik4) , dimension(klon) ::  ictop , icbot , ibotsc , idpl
    integer(ik4) , dimension(klon,klev) :: ilab
    logical , dimension(klon) :: ll_ldbase , llgo_on , lldeep , lldcum , &
             lldsc , llfirst , llresetjl
    logical :: llreset
    integer(ik4) :: icall , ik , ikb , is , jk , jl , jkk , jkt1 , jkt2 , jkt , jkb
    real(rk8) , dimension(klon,klev) :: zs , zsuh , zwu2h , zbuoh , zlu , &
             zqu , ztu , zuu , zvu , zcape
    real(rk8) , dimension(klon,klev+1) :: zsenh , zqenh
    real(rk8) ,dimension(klon) :: zqold , zph , zmix , zdz , zcbase , &
             ztven1 , ztvu1 , zdtvtrig
    real(rk8) :: zrho      ! density at surface (kg/m^3)
    real(rk8) :: zkhvfl    ! surface buoyancy flux (k m/s)
    real(rk8) :: zws       ! sigma_w at lowest model halflevel (m/s)
    real(rk8) :: zqexc     ! humidity excess at lowest model halflevel (kg/kg)
    real(rk8) :: ztexc     ! temperature excess at lowest model halflevel (k)
    real(rk8) :: zeps      ! fractional entrainment rate   [m^-1]
    real(rk8) :: ztvenh    ! environment virtual temperature at half levels (k)
    real(rk8) :: ztvuh     ! updraft virtual temperature at half levels     (k)
    real(rk8) :: zlglac    ! updraft liquid water frozen in one layer
    real(rk8) :: zqsu , zcor , zdq , zalfaw , zfacw , zfaci , zfac ,  &
                zesdp , zdqsdt , zdtdp , zdp ,zpdifftop ,zpdiffbot , &
                zsf , zqf , zbuof , zz , ztmp
    real(rk8) :: ztven2 , ztvu2 ! pseudoadiabatique t_v
    real(rk8) :: zwork1 , zwork2 ! work arrays for t and w perturbations

    real(rk8) , parameter :: zc2 = 0.55D0
    real(rk8) , parameter :: zaw = d_one
    real(rk8) , parameter :: zbw = d_one
    real(rk8) , parameter :: zepsadd = 1.D-4
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
                zwu2h(jl,jkk) = zws**2
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
                         min(d_one,(pqsen(jl,jk)/pqsen(jl,klev))**3)
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
        call cuadjtq(kidia,kfdia,klon,klev,ik,zph,ztu,zqu,llgo_on,icall)
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
              zfacw = c5les/((ztu(jl,ik)-c4les)**2)
              zfaci = c5ies/((ztu(jl,ik)-c4ies)**2)
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
            ! lldeep(jl) = paph(jl,jkb)-paph(jl,jkt) >= rdepths &
            !            .and. zdtvtrig(jl) > d_zero
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
!**** *CUCTRACER* - COMPUTE CONVECTIVE TRANSPORT OF CHEM. TRACERS
!                   IMPORTANT: ROUTINE IS FOR POSITIVE DEFINIT QUANTITIES
!
!          P.BECHTOLD        E.C.M.W.F.              11/02/2004
!
!**   INTERFACE.
!     ----------
!
!          *CUTRACER* IS CALLED FROM *CUMASTR*
!
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KTRAC*        NUMBER OF CHEMICAL TRACERS
!
!    *KTYPE*        CONVECTION TYPE (DEEP SHALLOW MID-LEVEL)
!    *KCTOP*        CLOUD TOP  LEVEL
!    *KCBOT*        CLOUD BASE LEVEL
!    *KDPL*         DEPARTURE LEVEL
!    *KDTOP*        DOWNDRAFT TOP LEVEL
!
!    INPUT PARAMETERS (LOGICAL):
!
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
!    *LDDRAF*       FLAG: .TRUE. IF DOWNDRAFTS EXIST
!
!    INPUT PARAMETERS (REAL):
!
!    *PTSPHY*       PHYSICS TIME-STEP                              S
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS            PA
!    *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATION
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PUDRATE       UPDRAFT DETRAINMENT                           KG/(M2*S)
!    *PDDRATE       DOWNDRAFT DETRAINMENT                         KG/(M2*S)
!
!    UPDATED PARAMETERS (REAL):
!
!    *PTENC*        UPDATED TENDENCY OF CHEM. TRACERS              1/S
!
!          METHOD
!          -------
!     EXPLICIT UPSTREAM AND IMPLICIT SOLUTION OF VERTICAL ADVECTION
!     DEPENDING ON VALUE OF RMFSOLCT: 0=EXPLICIT 0-1 SEMI-IMPLICIT >=1 IMPLICIT
!
!     FOR EXPLICIT SOLUTION: ONLY ONE SINGLE ITERATION
!     FOR IMPLICIT SOLUTION: FIRST IMPLICIT SOLVER, THEN EXPLICIT SOLVER
!                            TO CORRECT TENDENCIES BELOW CLOUD BASE
!
!
!------------------------------------------------------------------------------
!     COMMENTS FOR OFFLINE USERS IN CHEMICAL TRANSPORT MODELS
!     (i.e. reading mass fluxes and detrainment rates from ECMWF archive:
!      ------------------------------------------------------------------
!     KCTOP IS FIRST LEVEL FROM TOP WHERE PMFU>0
!     KDTOP IS FIRST LEVEL FROM TOP WHERE PMFD<0
!     KCBOT IS NOT NEEDED FOR EXPLICIT SOLUTION
!     ATTENTION: ON ARCHIVE DETRAINMENT RATES HAVE UNITS KG/(M3*S), SO FOR USE
!              IN CURRENT ROUTINE YOU HAVE TO MULTIPLY ARCHIVED VALUES BY DZ !!
!     LDCUM  IS TRUE IF CONVECTION EXISTS, i.e. IF PMFU>0 IN COLUMN OR IF
!                       KCTOP>0 AND KCTOP<KLEV
!     LDDRAF IS TRUE IF DOWNDRAUGHTS EXIST IF PMFD<0 IN COLUMN OR IF
!                       KDTOP>0 AND KDTOP<KLEV
!     IF MASSFLUX SATISFIES CFL CRITERIUM M<=DP/Dt IT IS SUFFICIENT TO
!     ONLY CONSIDER EXPLICIT SOLUTION (RMFSOLCT=0), IN THIS CASE
!     BUT IMPLICIT SOLUTION (RMFSOLCT=1) PART 7.0 IS PREFERED
!------------------------------------------------------------------------------
!
!          EXTERNALS
!          ---------
!          CUBIDIAG
!
!          MODIFICATIONS
!          -------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!
!----------------------------------------------------------------------
!
  subroutine cuctracer(kidia,kfdia,klon,ktdia,klev,ktrac,ktype,    &
                       kctop,kcbot,kdpl,kdtop,ldcum,lddraf,ptsphy, &
                       paph,pmfu,pmfd,pudrate,pddrate,pcen,ptenc)
    implicit none
    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: ktrac
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    integer(ik4) :: ktdia ! argument not used
    integer(ik4) , dimension(klon) , intent(in) :: ktype
    integer(ik4) , dimension(klon) , intent(in) :: kctop
    integer(ik4) , dimension(klon) , intent(in) :: kcbot
    integer(ik4) , dimension(klon) , intent(in) :: kdpl
    integer(ik4) , dimension(klon) , intent(in) :: kdtop
    logical , dimension(klon) , intent(in) :: ldcum
    logical , dimension(klon) , intent(in) :: lddraf
    real(rk8) , intent(in) :: ptsphy
    real(rk8) , dimension(klon,klev+1) , intent(in) :: paph
    real(rk8) , dimension(klon,klev) , intent(in) :: pmfu
    real(rk8) , dimension(klon,klev) , intent(in) :: pmfd
    real(rk8) , dimension(klon,klev) , intent(in) :: pudrate
    real(rk8) , dimension(klon,klev) , intent(in) :: pddrate
    real(rk8) , dimension(klon,klev,ktrac) , intent(in) :: pcen
    real(rk8) , dimension(klon,klev,ktrac) , intent(inout) :: ptenc

    !----------------------------------------------------------------------

    integer(ik4) :: ik , jk , jl , jn
    real(rk8) :: zzp , zmfa , zimp , zerate , zposi , ztsphy
    ! allocatable arays
    real(rk8) , dimension(:,:,:) , allocatable :: zcen , zcu , zcd , &
      ztenc , zmfc
    real(rk8) , dimension(:,:) , allocatable :: zdp , zb ,  zr1
    logical , dimension(:,:) ,  allocatable :: llcumask , llcumbas

    !----------------------------------------------------------------------

    zimp = d_one-rmfsolct
    ztsphy = d_one/ptsphy

    allocate(zcen(klon,klev,ktrac)) !half-level environmental values
    allocate(zcu(klon,klev,ktrac))  !updraft values
    allocate(zcd(klon,klev,ktrac))  !downdraft values
    allocate(ztenc(klon,klev,ktrac))!tendency
    allocate(zmfc(klon,klev,ktrac)) !fluxes
    allocate(zdp(klon,klev))        !pressure difference
    allocate(llcumask(klon,klev))   !mask for convection

    ! initialize cumulus mask + some setups

    do jk = 2 , klev
      do jl = kidia , kfdia
        llcumask(jl,jk) = .false.
        if ( ldcum(jl) ) then
          zdp(jl,jk) = egrav/(paph(jl,jk+1)-paph(jl,jk))
          if ( jk >= kctop(jl)-1 ) then
            llcumask(jl,jk) = .true.
          end if
        end if
      end do
    end do
    !----------------------------------------------------------------------

    do jn = 1 , ktrac

      !*    1.0          define tracers at half levels
      !                  -----------------------------
      do jk = 2 , klev
        ik = jk-1
        do jl = kidia , kfdia
          zcen(jl,jk,jn) = pcen(jl,jk,jn)
          zcd(jl,jk,jn) = pcen(jl,ik,jn)
          zcu(jl,jk,jn) = pcen(jl,ik,jn)
          zmfc(jl,jk,jn) = d_zero
          ztenc(jl,jk,jn) = d_zero
        end do
      end do
      do jl = kidia , kfdia
        zcu(jl,klev,jn) = pcen(jl,klev,jn)
      end do

      !*    2.0          compute updraft values
      !                  ----------------------
      do jk = klev-1 , 3 , -1
        ik = jk+1
        do jl = kidia , kfdia
          if ( llcumask(jl,jk) ) then
            zerate = pmfu(jl,jk)-pmfu(jl,ik)+pudrate(jl,jk)
            zmfa = d_one/max(rmfcmin,pmfu(jl,jk))
            if ( jk >= kctop(jl) )  then
              zcu(jl,jk,jn) = ( pmfu(jl,ik)*zcu(jl,ik,jn) + &
                zerate*pcen(jl,jk,jn)-pudrate(jl,jk)*zcu(jl,ik,jn) )*zmfa
              ! if you have a source term dc/dt=dcdt write
              ! zcu(jl,jk,jn) = ( pmfu(jl,ik)*zcu(jl,ik,jn) + &
              !   zerate*pcen(jl,jk,jn)-pudrate(jl,jk)*zcu(jl,ik,jn) )*zmfa + &
              !   dcdt(jl,ik,jn)*ptsphy
            end if
          end if
        end do
      end do

      !*    3.0          compute downdraft values
      !                  ------------------------
      do jk = 3 , klev
        ik = jk-1
        do jl = kidia , kfdia
          if ( lddraf(jl) .and. jk == kdtop(jl) ) then
            !nota: in order to avoid final negative tracer values at lfs
            ! the allowed value of zcd depends on the jump in mass flux
            ! at the lfs
            ! zcd(jl,jk,jn) = 0.5D0*zcu(jl,jk,jn)+0.5D0*pcen(jl,ik,jn)
            zcd(jl,jk,jn) = 0.1D0*zcu(jl,jk,jn)+0.9D0*pcen(jl,ik,jn)
          else if ( lddraf(jl) .and. jk > kdtop(jl) ) then
            zerate = -pmfd(jl,jk)+pmfd(jl,ik)+pddrate(jl,jk)
            zmfa = d_one/min(-rmfcmin,pmfd(jl,jk))
            zcd(jl,jk,jn) = ( pmfd(jl,ik)*zcd(jl,ik,jn) - &
              zerate*pcen(jl,ik,jn)+pddrate(jl,jk)*zcd(jl,ik,jn) )*zmfa
            ! if you have a source term dc/dt=dcdt write
            ! zcd(jl,jk,jn) = ( pmfd(jl,ik)*zcd(jl,ik,jn) - &
            !   zerate*pcen(jl,ik,jn)+pddrate(jl,jk)*zcd(jl,ik,jn) + &
            !   dcdt(jl,ik,jn)*ptsphy
          end if
        end do
      end do

      ! in order to avoid negative tracer at klev adjust zcd
      jk = klev
      ik = jk-1
      do jl = kidia , kfdia
        if ( lddraf(jl) ) then
          zposi = -zdp(jl,jk)*(pmfu(jl,jk)*zcu(jl,jk,jn) + &
                               pmfd(jl,jk)*zcd(jl,jk,jn) - &
                              (pmfu(jl,jk)+pmfd(jl,jk))*pcen(jl,ik,jn))
          if ( pcen(jl,jk,jn)+zposi*ptsphy < d_zero ) then
            zmfa = d_one/min(-rmfcmin,pmfd(jl,jk))
            zcd(jl,jk,jn) = ( (pmfu(jl,jk)+pmfd(jl,jk))*pcen(jl,ik,jn) - &
                               pmfu(jl,jk)*zcu(jl,jk,jn) + &
                               pcen(jl,jk,jn)/(ptsphy*zdp(jl,jk)) )*zmfa
          end if
        end if
      end do
    end do

    !----------------------------------------------------------------------

    do jn = 1 , ktrac

      !*    4.0          compute fluxes
      !                  --------------
      do jk = 2 , klev
        ik = jk-1
        do jl = kidia , kfdia
          if ( llcumask(jl,jk) ) then
            zmfa = pmfu(jl,jk)+pmfd(jl,jk)
            zmfc(jl,jk,jn) = pmfu(jl,jk)*zcu(jl,jk,jn) + &
                             pmfd(jl,jk)*zcd(jl,jk,jn) - &
                             zimp*zmfa*zcen(jl,ik,jn)
          end if
        end do
      end do

      !*    5.0          compute tendencies = rhs
      !                  ------------------------
      do jk = 2 , klev-1
        ik = jk+1
        do jl = kidia , kfdia
          if ( llcumask(jl,jk) ) then
            ztenc(jl,jk,jn) = zdp(jl,jk)*(zmfc(jl,ik,jn)-zmfc(jl,jk,jn))
          end if
        end do
      end do
      jk = klev
      do jl = kidia , kfdia
        if ( ldcum(jl) ) then
          ztenc(jl,jk,jn) = -zdp(jl,jk)*zmfc(jl,jk,jn)
        end if
      end do
    end do

    if ( rmfsolct == d_zero ) then
      !
      !*    6.0          update tendencies
      !                  -----------------
      do jn = 1 , ktrac
        do jk = 2 , klev
          do jl = kidia , kfdia
            if ( llcumask(jl,jk) ) then
              ptenc(jl,jk,jn) = ptenc(jl,jk,jn)+ztenc(jl,jk,jn)
            end if
          end do
        end do
      end do
    else
      !
      !*    7.0          implicit solution
      !                  -----------------
      ! fill bi-diagonal matrix vectors a=k-1, b=k;
      ! reuse zmfc=a and zb=b;
      ! ztenc corresponds to the rhs ("constants") of the equation
      ! the solution is in zr1
      allocate(zb(klon,klev))
      allocate(zr1(klon,klev))
      allocate(llcumbas(klon,klev))
      llcumbas(:,:) = .false.
      zb(:,:) = d_one

      do jn = 1 , ktrac
        ! fill vectors a, b and rhs
        do jk = 2 , klev
          ik = jk+1
          do jl = kidia , kfdia
            ! llcumbas(jl,jk) = llcumask(jl,jk).and.jk<=kcbot(jl)
            llcumbas(jl,jk) = llcumask(jl,jk)
            if ( llcumbas(jl,jk) ) then
              zzp = rmfsolct*zdp(jl,jk)*ptsphy
              zmfc(jl,jk,jn) = -zzp*(pmfu(jl,jk)+pmfd(jl,jk))
              ztenc(jl,jk,jn) = ztenc(jl,jk,jn)*ptsphy+pcen(jl,jk,jn)
              ! for implicit solution including tendency source term
              ! ztenc(jl,jk,jn) = (ztenc(jl,jk,jn) + &
              !       ptenc(jl,jk,jn))*ptsphy+pcen(jl,jk,jn)
              if ( jk < klev ) then
                zb(jl,jk) = d_one+zzp*(pmfu(jl,ik)+pmfd(jl,ik))
              else
                zb(jl,jk) = d_one
              end if
            end if
          end do
        end do
        call cubidiag(kidia,kfdia,klon,klev,kctop,llcumbas, &
                      zmfc(:,:,jn),zb,ztenc(:,:,jn),zr1)
        ! compute tendencies
        do jk = 2 , klev
          do jl = kidia , kfdia
            if ( llcumbas(jl,jk) ) then
              ptenc(jl,jk,jn) = ptenc(jl,jk,jn) + &
                (zr1(jl,jk)-pcen(jl,jk,jn))*ztsphy
              !  for implicit solution including tendency source term
              !  ptenc(jl,jk,jn) = (zr1(jl,jk)-pcen(jl,jk,jn))*ztsphy
            end if
          end do
        end do
      end do
      deallocate(llcumbas)
      deallocate(zb)
      deallocate(zr1)
    end if

    !---------------------------------------------------------------------------

    deallocate(llcumask)
    deallocate(zdp)
    deallocate(zmfc)
    deallocate(ztenc)
    deallocate(zcd)
    deallocate(zcu)
    deallocate(zcen)
  end subroutine cuctracer
!
!**** *CUMASTR*  MASTER ROUTINE FOR CUMULUS MASSFLUX-SCHEME
!
!     M.TIEDTKE      E.C.M.W.F.     1986/1987/1989
!     D.GREGORY      E.C.M.W.F.     1996
!     P.BECHTOLD     E.C.M.W.F.     2005/2007
!
!     PURPOSE
!     -------
!
!          THIS ROUTINE COMPUTES THE PHYSICAL TENDENCIES OF THE
!     PROGNOSTIC VARIABLES T,Q,U, V AND TRACERS DUE TO CONVECTIVE PROCESSES.
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
!      3.CLOUD BASE, CLOUD TOP AND PRECIP FOR RADIATION
!        (USED IN SUBROUTINE CLOUD)
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
!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):
!
!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KTDIA*        START OF THE VERTICAL LOOP
!    *KLEV*         NUMBER OF LEVELS
!    *KTRAC*        NUMBER OF CHEMICAL TRACERS
!    *KSTEP*        CURRENT TIME STEP INDEX
!    *KSTART*       FIRST STEP OF MODEL
!
!     INPUT PARAMETERS (LOGICAL)
!
!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)
!
!     INPUT PARAMETERS (REAL)
!
!    *PTSPHY*       TIME STEP FOR THE PHYSICS                       S
!    *PTEN*         PROVISIONAL ENVIRONMENT TEMPERATURE (T+1)       K
!    *PQEN*         PROVISIONAL ENVIRONMENT SPEC. HUMIDITY (T+1)  KG/KG
!    *PUEN*         PROVISIONAL ENVIRONMENT U-VELOCITY (T+1)       M/S
!    *PVEN*         PROVISIONAL ENVIRONMENT V-VELOCITY (T+1)       M/S
!    *PCEN*         PROVISIONAL ENVIRONMENT TRACER CONCENTRATIONS KG/KG
!    *PLITOT*       GRID MEAN LIQUID WATER+ICE CONTENT            KG/KG
!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
!    *PQSEN*        ENVIRONMENT SPEC. SATURATION HUMIDITY (T+1)   KG/KG
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)        KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS             PA
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS             PA
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!
!    UPDATED PARAMETERS (REAL):
!
!    *PTENT*        TEMPERATURE TENDENCY                           K/S
!    *PTENQ*        MOISTURE TENDENCY                             KG/(KG S)
!    *PTENL*        LIQUID WATER TENDENCY                         KG/(KG S)
!    *PTENI*        ICE CONDENSATE TENDENCY                       KG/(KG S)
!    *PTENU*        TENDENCY OF U-COMP. OF WIND                    M/S2
!    *PTENV*        TENDENCY OF V-COMP. OF WIND                    M/S2
!    *PTENC*        TENDENCY OF CHEMICAL TRACERS                   1/S
!
!    OUTPUT PARAMETERS (LOGICAL):
!
!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS
!    *LDSC*         FLAG: .TRUE. FOR SC-POINTS
!
!    OUTPUT PARAMETERS (INTEGER):
!
!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL
!    *KBOTSC*       CLOUD BASE LEVEL FOR SC-CLOUDS
!
!    OUTPUT PARAMETERS (REAL):
!
!    *PTU*          TEMPERATURE IN UPDRAFTS                         K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PLUDE*        DETRAINED LIQUID WATER                        KG/(M2*S)
!    *PENTH*        INCREMENT OF DRY STATIC ENERGY                 J/(KG*S)
!    *PMFLXR*       CONVECTIVE RAIN FLUX                          KG/(M2*S)
!    *PMFLXS*       CONVECTIVE SNOW FLUX                          KG/(M2*S)
!    *PRAIN*        TOTAL PRECIP. PRODUCED IN CONV. UPDRAFTS      KG/(M2*S)
!                   (NO EVAPORATION IN DOWNDRAFTS)
!    *PMFU*         MASSFLUX UPDRAFTS                             KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                           KG/(M2*S)
!    *PMFUDE_RATE*  UPDRAFT DETRAINMENT RATE                      KG/(M3*S)
!    *PMFDDE_RATE*  DOWNDRAFT DETRAINMENT RATE                    KG/(M3*S)
!    *PCAPE*        CONVECTVE AVAILABLE POTENTIAL ENERGY           J/KG
!    *PWMEAN*       VERTICALLY AVERAGED UPDRAUGHT VELOCITY         M/S
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
!       CUCTRACER: TRACER TRANSPORT
!
!     SWITCHES.
!     --------
!
!          LMFPEN=.TRUE.   PENETRATIVE CONVECTION IS SWITCHED ON
!          LMFSCV=.TRUE.   SHALLOW CONVECTION IS SWITCHED ON
!          LMFMID=.TRUE.   MIDLEVEL CONVECTION IS SWITCHED ON
!          LMFDD=.TRUE.    CUMULUS DOWNDRAFTS SWITCHED ON
!          LMFDUDV=.TRUE.  CUMULUS FRICTION SWITCHED ON
!          LMFTRAC=.false. TRACER TRANSPORT
!
!
!     REFERENCE.
!     ----------
!
!          PAPER ON MASSFLUX SCHEME (TIEDTKE,1989)
!          DRAFT PAPER ON MASSFLUX SCHEME (NORDENG, 1995)
!          Bechtold et al. (2008 QJRMS 134,1337-1351), Rooy et al. (2012 QJRMS)
!
!          MODIFICATIONS
!          -------------
!             92-09-21 : Update to Cy44      J.-J. MORCRETTE
!             96-03-12 : Introduce CAPE closure for deep convection
!                        (based upon previous work by Nordeng, 1995)
!             99-06-21 : Optimisation   D.Salmond
!             03-08-29 : Clean-up, deep/shallow switches  P.Bechtold
!             04-02-11 : Add tracer transport             P.Bechtold
!             05-02-11 : Positive scaling of total Mflux  P.Bechtold
!             M.Hamrud      01-Oct-2003 CY28 Cleaning
!             04-12-03 : Turn off shallow convection over stratocu.
!                                                         M.Ko"hler
!             05-06-27 : Switch off ddraught if idtop<kctop
!                        correction for detrainment rates P.Bechtold
!             05-11-22 : Mods for coarser/finer physics D.Salmond + M.Hortal
!             06-02-11 : Enable TQ implicit               P.Bechtold
!             07-06-01 : Only single updraught call with  P.Bechtold
!                        scaling, convective turnover time
!                        scale, contain momentum computations
!                        in cumastrn
!             07-10-09 : Added KE dissipation             P. Bechtold
!             02-03-12 : remove all entrainment stuff   P. Bechtold
!----------------------------------------------------------------------
!
  subroutine cumastrn(kidia,kfdia,klon,ktdia,klev,kstep,kstart, &
                      ldland,ptsphy,pten,pqen,puen,pven,plitot, &
                      pvervel,pqhfl,pahfs,pap,paph,pgeo,pgeoh,  &
                      ptent,ptenq,ptenu,ptenv,ptenl,pteni,      &
                      ldcum,ktype,kcbot,kctop,kbotsc,ldsc,      &
                      ptu,pqu,plu,pmflxr,pmflxs,prain,pmfu,pmfd,&
                      pmfude_rate,pmfdde_rate,pcape,ktrac,pcen, &
                      ptenc)
    implicit none

    integer(ik4) , intent(in) :: klon
    integer(ik4) , intent(in) :: klev
    integer(ik4) , intent(in) :: kidia
    integer(ik4) , intent(in) :: kfdia
    integer(ik4) , intent(in) :: ktrac
    integer(ik4) :: ktdia ! argument not used
    integer(ik4) :: kstep ! argument not used
    integer(ik4) :: kstart ! argument not used
    logical , dimension(klon) , intent(in) :: ldland
    real(rk8) , intent(in)    :: ptsphy
    real(rk8) , dimension(klon,klev) , intent(inout) :: pten
    real(rk8) , dimension(klon,klev) , intent(inout) :: pqen
    real(rk8) , dimension(klon,klev) , intent(in) :: puen
    real(rk8) , dimension(klon,klev) , intent(in) :: pven
    real(rk8) , dimension(klon,klev) , intent(in) :: plitot
    real(rk8) , dimension(klon,klev) , intent(in) :: pvervel
    !real(rk8) , dimension(klon,klev) , intent(inout) :: pqsen
    real(rk8) , dimension(klon,klev+1) , intent(in) :: pqhfl
    real(rk8) , dimension(klon,klev+1) , intent(in) :: pahfs
    real(rk8) , dimension(klon,klev) , intent(in) :: pap
    real(rk8) , dimension(klon,klev+1) , intent(in) :: paph
    real(rk8) , dimension(klon,klev) , intent(in) :: pgeo
    real(rk8) , dimension(klon,klev+1) , intent(in) :: pgeoh
    real(rk8) , dimension(klon,klev,ktrac) , intent(in) :: pcen
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptent
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptenq
    real(rk8) , dimension(klon,klev) , intent(out)   :: ptenl
    real(rk8) , dimension(klon,klev) , intent(out)   :: pteni
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptenu
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptenv
    real(rk8) , dimension(klon,klev,ktrac) , intent(inout) :: ptenc
    logical , dimension(klon) , intent(inout) :: ldcum
    integer(ik4) , dimension(klon) , intent(inout) :: ktype
    integer(ik4) , dimension(klon) , intent(inout) :: kcbot
    integer(ik4) , dimension(klon) , intent(inout) :: kctop
    integer(ik4) , dimension(klon) , intent(out)   :: kbotsc
    logical , dimension(klon) , intent(out)   :: ldsc
    !logical , dimension(klon) , intent(in)    :: ldshcv
    real(rk8) , dimension(klon,klev) , intent(inout) :: ptu
    real(rk8) , dimension(klon,klev) , intent(inout) :: pqu
    real(rk8) , dimension(klon,klev) , intent(inout) :: plu
    !real(rk8) , dimension(klon,klev) , intent(inout) :: plude
    !real(rk8) , dimension(klon,klev) , intent(out)   :: penth
    real(rk8) , dimension(klon,klev+1) , intent(inout) :: pmflxr
    real(rk8) , dimension(klon,klev+1) , intent(inout) :: pmflxs
    real(rk8) , dimension(klon) , intent(out) :: prain
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfu
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfd
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfude_rate
    real(rk8) , dimension(klon,klev) , intent(inout) :: pmfdde_rate
    real(rk8) , dimension(klon) , intent(out) :: pcape
    !real(rk8) , intent(out)   :: pwmean(klon)
    !*upg change to operations
    real(rk8) , dimension(klon) :: pwmean
    real(rk8) , dimension(klon,klev) :: plude ! only local variable
    real(rk8) , dimension(klon,klev) :: penth ! only local variable
    real(rk8) , dimension(klon,klev) :: pqsen ! only local variable
    !*upg
    real(rk8) , dimension(klon,klev) :: ztenh , zqenh , zqsenh , &
            ztd , zqd , zmfus , zmfds , zmfuq , zmfdq , zdmfup , &
            zdmfdp , zmful , zuu , zvu , zud , zvd , zkineu , zkined
    real(rk8) , dimension(klon) :: zrfl
    real(rk8) , dimension(klon) :: zmfub , zmfub1 , zdqcv
    real(rk8) , dimension(klon,klev) :: zdpmel , zlglac
    real(rk8) , dimension(klon) :: zdhpbl , zwubase
    real(rk8) , dimension(klon,klev) :: zdmfen , zdmfde
    integer(ik4) , dimension(klon,klev) ::  ilab
    integer(ik4) , dimension(klon) :: idtop , ictop0 , ilwmin
    integer(ik4) , dimension(klon) :: idpl ! departure level for convection
    real(rk8) , dimension(klon) :: zcape , zheat
    logical , dimension(klon) :: llddraf , llddraf3 , lldcum
    logical ::  llo1
    logical , dimension(klon) :: llo2
    integer(ik4) :: ikb, itopm2, jk, ik, jl
    real(rk8) ::   zcons2 , zcons , zdh , zdqmin , zdz , zeps , zfac , &
            zmfmax , zpbmpt , zqumqe , zro , zmfa , zerate , zderate , &
            zduten , zdvten , ztdis , zalv
    real(rk8) , dimension(klon) :: zsfl
    real(rk8) , dimension(klon) :: ztau ! adjustment time

    ! scaling factor for momentum and tracer massflux
    real(rk8) , dimension(klon) :: zmfs , zmfuub , zmfuvb , zsum12 , &
      zsum22 , zmf_shal
    real(rk8) , dimension(klon,klev) :: zmfuus , zmfdus , zmfudr , &
      zmfddr , ztenu , ztenv , zuv2

    !*upg change to operations
    !    locals for conservation check
    logical :: llconscheck = .false.
    integer(ik4) :: jn
    real(rk8) , dimension(:,:) , allocatable :: ztent , ztenq , zsumc
    real(rk8) , dimension(:,:,:) , allocatable :: ztenc

    !---------------------------------------------------------------------
    !*upg change to operations call satur routine here
    !     0.           compute saturation specific humidity
    !                  ------------------------------------
    ldcum(:) = .false.
    pqsen(:,:) = pqen(:,:)
    call satur(kidia,kfdia,klon,njkt2,klev,pap,pten,pqsen,1)

    !*upg
    !---------------------------------------------------------------------
    !     1.           specify constants and parameters
    !                  --------------------------------

    zcons2 = rmfcfl/(egrav*ptsphy)
    zcons = d_one/(egrav*ptsphy)

    !*    2.           initialize values at vertical grid points in 'cuini'
    !                  ---------------------------------------------------

    call cuinin(kidia,kfdia,klon,klev,pten,pqen,pqsen,puen,pven,  &
                pvervel,pgeo,paph,ilwmin,ilab,ztenh,zqenh,zqsenh, &
                pgeoh,ptu,pqu,ztd,zqd,zuu,zvu,zud,zvd,plu)

    !---------------------------------------------------------------------
    !*    3.0          cloud base calculations
    !                  -----------------------
    !*             (a) determine cloud base values in 'cubase'
    !                  ---------------------------------------

    call cubasen(kidia,kfdia,klon,ktdia,klev,ldland,ztenh,zqenh,pgeoh, &
                 paph,pqhfl,pahfs,pten,pqen,pqsen,pgeo,puen,pven,      &
                 ptu,pqu,plu,zuu,zvu,zwubase,ilab,ldcum,ldsc,kcbot,    &
                 kbotsc,ictop0,idpl,pcape)

    !*             (b) determine total moisture convergence and
    !*                 decide on type of cumulus convection
    !*                 one the basis of the depth of the convection
    !*                 deep if cloud depth > 200mb
    !*                 shallow if cloud depth <200mb
    !                  -----------------------------------------
    ! calculate column and sub cloud layer moisture convergence
    ! and sub cloud layer moist static energy convergence

    do jl = kidia , kfdia
      zdqcv(jl)  = d_zero
      zdhpbl(jl) = d_zero
      idtop(jl) = 0
    end do
    do jk = njkt2 , klev
      do jl = kidia , kfdia
        zdqcv(jl)=zdqcv(jl)+max(d_zero,ptenq(jl,jk))*(paph(jl,jk+1)-paph(jl,jk))
        if ( ldcum(jl) .and. jk >= kcbot(jl) ) then
          zdhpbl(jl)=zdhpbl(jl)+(wlhv*ptenq(jl,jk)+cpd*ptent(jl,jk)) * &
                     (paph(jl,jk+1)-paph(jl,jk))
        end if
      end do
    end do

    !*                 estimate cloud height for entrainment/detrainment
    !*                 calculations in cuasc and initial determination of
    !*                 cloud type
    !*                 (max.possible cloud height
    !*                 for non-entraining plume, following a.-s.,1974)
    !                  -------------------------------------------------

!dir$ ivdep
!ocl novrec

    !*                 specify initial cloud type
    !*

    do jl = kidia , kfdia
      if ( ldcum(jl) ) then
        ikb = kcbot(jl)
        itopm2 = ictop0(jl)
        zpbmpt = paph(jl,ikb)-paph(jl,itopm2)
        if ( zpbmpt >= rdepths ) then
          ktype(jl) = 1
        else
          ktype(jl) = 2
        end if
      else
        ktype(jl) = 0
      end if
    end do

    !*             (c) calculate initial updraught mass flux
    !*                 and set lateral mixing rates
    !*
    !*                 for deep convection assume it is 10% of
    !*                 maximum value which is determined by the
    !*                 thickness of the layer and timestep
    !*
    !*                 for shallow convection calculated assuming
    !*                 a balance of moist static energy in the
    !*                 sub-cloud layer (ignores present of downdraughts)
    !                  ------------------------------------------

!dir$ ivdep
!ocl novrec
    if ( lmfwstar ) then
      do jl = kidia , kfdia
        if ( ldcum(jl) ) then
          ikb = kcbot(jl)
          zdz = max(d_zero,min(1.5D3,(pgeoh(jl,ikb)-pgeoh(jl,klev+1))/egrav))
          zmf_shal(jl) = 0.07D0*(egrav/pten(jl,klev)*zdz * &
                         max(d_zero,-pahfs(jl,klev+1)*rcpd-retv * &
                             pten(jl,klev)*pqhfl(jl,klev+1)))**0.3333D0
          zmfmax = (paph(jl,ikb)-paph(jl,ikb-1))*zcons2*rmflic+rmflia
          zmf_shal(jl) = min(zmf_shal(jl),zmfmax)
        end if
      end do
    end if

    do jl = kidia , kfdia
      if ( ldcum(jl) ) then
        ikb = kcbot(jl)
        zmfmax = (paph(jl,ikb)-paph(jl,ikb-1))*zcons2*rmflic+rmflia

        ! deep convection

        if ( ktype(jl) == 1 ) then
          zmfub(jl) = zmfmax*0.1D0

        else if ( ktype(jl) == 2 ) then

        ! shallow convection

          zqumqe = pqu(jl,ikb)+plu(jl,ikb)-zqenh(jl,ikb)
          zdqmin = max(0.01D0*zqenh(jl,ikb),1.0D-10)
          zdh = cpd*(ptu(jl,ikb)-ztenh(jl,ikb))+wlhv*zqumqe
          zdh = egrav*max(zdh,1.D5*zdqmin)
          if ( zdhpbl(jl) > d_zero ) then
            zmfub(jl) = zdhpbl(jl)/zdh
            !eps: temporary solution for explicit
            if( ptsphy > 1800.D0 .and. rmfcfl == d_one ) then
              zmfub(jl) = min(zmfub(jl),3.D0*zmfmax)
            else
              zmfub(jl) = min(zmfub(jl),zmfmax)
            end if
          else
            zmfub(jl) = zmfmax*0.1D0
            ldcum(jl) = .false.
          end if
          if ( lmfwstar ) zmfub(jl) = zmf_shal(jl)
        end if

      else

        ! no buoyancy cloud base from surface
        ! set cloud base mass flux and mixing rate
        ! to default value for safety

        zmfub(jl) = d_zero
      end if
    end do

    !-----------------------------------------------------------------------
    !*    4.0          determine cloud ascent for entraining plume
    !                  -------------------------------------------
    !*             (a) estimate cloud height for entrainment/detrainment
    !*                 calculations in cuasc (max.possible cloud height
    !*                 for non-entraining plume, following a.-s.,1974)
    !                  -------------------------------------------------
    ! calculations now done is section 3 above so that
    ! initial cloud depth can be used to specify
    ! the type of convection

    !*             (b) do ascent in 'cuasc'in absence of downdrafts
    !                  --------------------------------------------

    call cuascn(kidia,kfdia,klon,klev,ptsphy,ztenh,zqenh,pten,  &
                pqen,pqsen,plitot,pgeo,pgeoh,pap,paph,pvervel,  &
                zwubase,ldland,ldcum,ktype,ilab,ptu,pqu,plu,    &
                pmfu,zmfub,zlglac,zmfus,zmfuq,zmful,plude,      &
                zdmfup,zdmfen,kcbot,kctop,ictop0,idpl,          &
                pmfude_rate,zkineu,pwmean)

    !*         (c) check cloud depth and change entrainment rate accordingly
    !              calculate precipitation rate (for downdraft calculation)
    !              -----------------------------------------------------
!dir$ ivdep
!ocl novrec
    do jl = kidia , kfdia
      if ( ldcum(jl) ) then
        ikb = kcbot(jl)
        itopm2 = kctop(jl)
        zpbmpt = paph(jl,ikb)-paph(jl,itopm2)
        if ( ktype(jl) == 1 .and. zpbmpt < rdepths ) ktype(jl) = 2
        if ( ktype(jl) == 2 .and. zpbmpt >= rdepths ) ktype(jl) = 1
        ictop0(jl) = kctop(jl)
      end if
      zrfl(jl) = zdmfup(jl,1)
    end do
    do jk = 2 , klev
      do jl = kidia , kfdia
        zrfl(jl) = zrfl(jl)+zdmfup(jl,jk)
      end do
    end do
    do jk = 1 , klev
      do jl = kidia , kfdia
        pmfd(jl,jk) = d_zero
        zmfds(jl,jk) = d_zero
        zmfdq(jl,jk) = d_zero
        zdmfdp(jl,jk) = d_zero
        zdpmel(jl,jk) = d_zero
      end do
    end do
    if ( lmfuvdis ) then
      do jk = 1 , klev
        do jl = kidia , kfdia
          ztenu(jl,jk) = ptenu(jl,jk)
          ztenv(jl,jk) = ptenv(jl,jk)
        end do
      end do
    end if

    !-----------------------------------------------------------------------
    !*    5.0          cumulus downdraft calculations
    !                  ------------------------------

    if ( lmfdd ) then

      !*             (a) determine lfs in 'cudlfs'
      !                  -------------------------

      call cudlfsn(kidia,kfdia,klon,klev,kcbot,kctop,ldcum, &
                   ztenh,zqenh,pten,pqsen,pgeo,pgeoh,paph,  &
                   ptu,pqu,zmfub,zrfl,ztd,zqd,pmfd,zmfds,   &
                   zmfdq,zdmfdp,idtop,llddraf)

      !*            (b)  determine downdraft t,q and fluxes in 'cuddraf'
      !                  -----------------------------------------------

      call cuddrafn(kidia,kfdia,klon,klev,llddraf,ztenh,zqenh,    &
                    pgeo,pgeoh,paph,zrfl,ztd,zqd,pmfu,pmfd,zmfds, &
                    zmfdq,zdmfdp,zdmfde,pmfdde_rate,zkined)

    end if

    !-----------------------------------------------------------------------
    !*    6.0          closure
    !                  ------
    !*                 recalculate cloud base massflux from a
    !*                 cape closure for deep convection (ktype=1)
    !*                 and by pbl equilibrum taking downdrafts into
    !*                 account for shallow convection (ktype=2)
    !                  --------------------------------------------

    !   deep convection

    do jl = kidia , kfdia
      zheat(jl) = d_zero
      zcape(jl) = d_zero
      zmfub1(jl) = zmfub(jl)
    end do
    do jk = 1 , klev
      do jl = kidia , kfdia
        llo1 = ldcum(jl) .and. ktype(jl) == 1
        if ( llo1 .and. jk <= kcbot(jl) .and. jk > kctop(jl) ) then
          ikb = kcbot(jl)
          zro = paph(jl,jk)/(rgas*ztenh(jl,jk)*(d_one+retv*zqenh(jl,jk)))
          zdz = (pgeoh(jl,jk-1)-pgeoh(jl,jk))
          zheat(jl) = zheat(jl) + &
                (  (pten(jl,jk-1)-pten(jl,jk) + zdz*rcpd)/ztenh(jl,jk) + &
                    retv*(pqen(jl,jk-1)-pqen(jl,jk))  ) * &
                (egrav*(pmfu(jl,jk)+pmfd(jl,jk)))/zro
          zcape(jl) = zcape(jl) + &
                 ((ptu(jl,jk)-ztenh(jl,jk))/ztenh(jl,jk) + &
                   retv*(pqu(jl,jk)-zqenh(jl,jk))-plu(jl,jk) ) * zdz
        end if
      end do
    end do

    do jl = kidia , kfdia
      if ( ldcum(jl) .and. ktype(jl) == 1 ) then
        ikb = kcbot(jl)
        ik = kctop(jl)
        zcape(jl) = max(d_zero,min(zcape(jl),5000.0D0))
        zheat(jl) = max(1.D-4,zheat(jl))
        ztau(jl) = (pgeoh(jl,ik)-pgeoh(jl,ikb)) / &
          ((2.0D0+min(15.0D0,pwmean(jl)))*egrav)*rtau
        ztau(jl) = max(ptsphy,min(10800.D0,ztau(jl)))
        ztau(jl) = max(720.0D0,ztau(jl))
        zmfub1(jl) = (zcape(jl)*zmfub(jl))/(zheat(jl)*ztau(jl))
        zmfub1(jl) = max(zmfub1(jl),0.001D0)
        zmfmax = (paph(jl,ikb)-paph(jl,ikb-1))*zcons2*rmflic+rmflia
        zmfub1(jl) = min(zmfub1(jl),zmfmax)
      end if
    end do

    !  shallow convection and mid_level

!dir$ ivdep
!ocl novrec
    do jl = kidia , kfdia
      if ( ldcum(jl) .and. (ktype(jl) == 2 .or. ktype(jl) == 3) ) then
        ikb = kcbot(jl)
        if ( pmfd(jl,ikb) < d_zero ) then
          zeps = -pmfd(jl,ikb)/max(zmfub(jl),1.D-10)
        else
          zeps = d_zero
        end if
        zqumqe = pqu(jl,ikb)+plu(jl,ikb) - &
                 zeps*zqd(jl,ikb)-(d_one-zeps)*zqenh(jl,ikb)
        zdqmin = max(0.01D0*zqenh(jl,ikb),1.D-10)
        ! maximum permisable value of ud base mass flux
        zmfmax = (paph(jl,ikb)-paph(jl,ikb-1))*zcons2*rmflic+rmflia

        ! shallow convection

        if ( ktype(jl) == 2 ) then
          zdh = cpd*(ptu(jl,ikb)-zeps*ztd(jl,ikb) - &
               (d_one-zeps)*ztenh(jl,ikb))+wlhv*zqumqe
          zdh = egrav*max(zdh,1.D5*zdqmin)
          if ( zdhpbl(jl) > d_zero ) then
            zmfub1(jl) = zdhpbl(jl)/zdh
          else
            zmfub1(jl) = zmfub(jl)
          end if
          !eps: temporary solution for explicit
          if ( ptsphy>1800.D0 .and. rmfcfl==d_one ) then
            zmfub1(jl) = min(zmfub1(jl),3.D0*zmfmax)
          else
            zmfub1(jl) = min(zmfub1(jl),zmfmax)
          end if
          if ( lmfwstar ) zmfub1(jl) = zmf_shal(jl)
        end if

        ! mid-level convection

        if ( ktype(jl) == 3 ) then
          zmfub1(jl) = zmfub(jl)*(d_one+zeps)
          zmfub1(jl) = min(zmfub1(jl),zmfmax)
        end if

      end if
    end do

    ! rescale dd fluxes if deep and shallow convection

    do jk = 1 , klev
      do jl = kidia , kfdia
        if ( llddraf(jl) .and. ( ktype(jl) == 1 .or. ktype(jl) == 2 ) ) then
          zfac = zmfub1(jl)/max(zmfub(jl),1.D-10)
          pmfd(jl,jk) = pmfd(jl,jk)*zfac
          zmfds(jl,jk) = zmfds(jl,jk)*zfac
          zmfdq(jl,jk) = zmfdq(jl,jk)*zfac
          zdmfdp(jl,jk) = zdmfdp(jl,jk)*zfac
          !  also rescale detrainment flux for era pp
          pmfdde_rate(jl,jk) = pmfdde_rate(jl,jk)*zfac
        end if
      end do
    end do

    !-----------------------------------------------------------------------
    !*    6.2          final closure=scaling
    !                  ---------------------

    do jl = kidia , kfdia
      if ( ldcum(jl) ) then
        zmfs(jl) = zmfub1(jl)/max(rmfcmin,zmfub(jl))
      end if
    end do
    do jk = 2 , klev
      do jl = kidia , kfdia
        if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
          ikb = kcbot(jl)
          if ( jk>ikb ) then
            zdz = ((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ikb)))
            pmfu(jl,jk) = pmfu(jl,ikb)*zdz
          end if
          zmfmax = (paph(jl,jk)-paph(jl,jk-1))*zcons2*rmflic+rmflia
          if ( pmfu(jl,jk)*zmfs(jl)>zmfmax ) then
            zmfs(jl) = min(zmfs(jl),zmfmax/pmfu(jl,jk))
          end if
        end if
      end do
    end do
    do jk = 2 , klev
      do jl = kidia , kfdia
        if ( ldcum(jl) .and. jk <= kcbot(jl) .and. jk >= kctop(jl)-1 ) then
          pmfu(jl,jk) = pmfu(jl,jk)*zmfs(jl)
          zmfus(jl,jk) = zmfus(jl,jk)*zmfs(jl)
          zmfuq(jl,jk) = zmfuq(jl,jk)*zmfs(jl)
          zmful(jl,jk) = zmful(jl,jk)*zmfs(jl)
          zdmfup(jl,jk) = zdmfup(jl,jk)*zmfs(jl)
          zdmfen(jl,jk) = zdmfen(jl,jk)*zmfs(jl)
          plude(jl,jk) = plude(jl,jk)*zmfs(jl)
          pmfude_rate(jl,jk) = pmfude_rate(jl,jk)*zmfs(jl)
        end if
      end do
    end do

    !-----------------------------------------------------------------------
    !*    6.5          in case that either deep or shallow is switched off
    !                  reset ldcum to false-> fluxes set to zero in cuflxn
    !                  ---------------------------------------------------
    !                 exclude pathological ktype=2 kcbot=kctop=klev-1

    do jl = kidia , kfdia
      if ( ktype(jl) == 2 .and. &
           kcbot(jl) == kctop(jl) .and. kcbot(jl) >= klev-1 ) then
        ldcum(jl) = .false.
        ktype(jl) = 0
      end if
    end do

    if ( .not. lmfscv .or. .not. lmfpen ) then
      do jl = kidia , kfdia
        llo2(jl) = .false.
        if ( ( .not. lmfscv .and. ktype(jl) == 2 ) .or. &
             ( .not. lmfpen .and. ktype(jl) == 1) )then
          llo2(jl) = .true.
          ldcum(jl) = .false.
        end if
      end do
    end if

    !-----------------------------------------------------------------------
    !*    7.0          determine final convective fluxes in 'cuflx'
    !                  ------------------------------------------
    !- set dd mass fluxes to zero above cloud top
    !  (because of inconsistency with second updraught)
    do jl = kidia , kfdia
      if ( llddraf(jl) .and. idtop(jl) <= kctop(jl) ) then
        idtop(jl) = kctop(jl)+1
      end if
    end do

    do jk = 2 , klev
      do jl = kidia , kfdia
        if ( llddraf(jl) ) then
          if ( jk < idtop(jl) ) then
            pmfd(jl,jk) = d_zero
            zmfds(jl,jk) = d_zero
            zmfdq(jl,jk) = d_zero
            pmfdde_rate(jl,jk) = d_zero
            zdmfdp(jl,jk) = d_zero
          else if ( jk == idtop(jl) ) then
            pmfdde_rate(jl,jk) = d_zero
          end if
        end if
      end do
    end do
    call cuflxn(kidia,kfdia,klon,klev,ptsphy,pten,pqen,pqsen,ztenh,  &
                zqenh,paph,pap,pgeoh,ldland,ldcum,kcbot,kctop,idtop, &
                itopm2,ktype,llddraf,pmfu,pmfd,zmfus,zmfds,zmfuq,    &
                zmfdq,zmful,plude,zdmfup,zdmfdp,zdpmel,zlglac,pmflxr,&
                pmflxs,prain,pmfdde_rate)

    !- rescale dd fluxes if total mass flux becomes negative
    !- correct dd detrainment rates if entrainment becomes negative
    !- correct ud detrainment rates if entrainment becomes negative
    !- conservation correction for precip

    zmfs(:) = d_one
    ! do jk = 2 , klev-1
    do jk = 2 , klev ! change for stability
      do jl = kidia , kfdia
        if ( llddraf(jl) .and. jk >= idtop(jl)-1 ) then
          zmfmax = pmfu(jl,jk)*0.98D0
          if ( pmfd(jl,jk)+zmfmax+1.D-15 < d_zero ) then
            zmfs(jl) = min(zmfs(jl),-zmfmax/pmfd(jl,jk))
          end if
        end if
      end do
    end do
    zmfuub(:) = d_zero
    do jk = 2 , klev
      do jl = kidia , kfdia
        if ( zmfs(jl) < d_one .and. jk >= idtop(jl)-1 ) then
          pmfd(jl,jk) = pmfd(jl,jk)*zmfs(jl)
          zmfds(jl,jk) = zmfds(jl,jk)*zmfs(jl)
          zmfdq(jl,jk) = zmfdq(jl,jk)*zmfs(jl)
          pmfdde_rate(jl,jk) = pmfdde_rate(jl,jk)*zmfs(jl)
          zmfuub(jl) = zmfuub(jl)-(d_one-zmfs(jl))*zdmfdp(jl,jk)
          pmflxr(jl,jk+1) = pmflxr(jl,jk+1)+zmfuub(jl)
          zdmfdp(jl,jk) = zdmfdp(jl,jk)*zmfs(jl)
        end if
      end do
    end do

    do jk = 2 , klev-1
      do jl = kidia , kfdia
        if ( llddraf(jl) .and. jk >= idtop(jl)-1 ) then
          zerate = -pmfd(jl,jk)+pmfd(jl,jk-1)+pmfdde_rate(jl,jk)
          if ( zerate < d_zero ) then
            pmfdde_rate(jl,jk) = pmfdde_rate(jl,jk)-zerate
          end if
        end if
        if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
          zerate = pmfu(jl,jk)-pmfu(jl,jk+1)+pmfude_rate(jl,jk)
          if ( zerate < d_zero ) then
            pmfude_rate(jl,jk) = pmfude_rate(jl,jk)-zerate
          end if
          ! zdmfup(jl,jk) = zdmfup(jl,jk)+zdmfdp(jl,jk)
          zdmfup(jl,jk) = pmflxr(jl,jk+1)+pmflxs(jl,jk+1) - &
                          pmflxr(jl,jk)-pmflxs(jl,jk)
          zdmfdp(jl,jk) = d_zero
        end if
      end do
    end do

    ! avoid negative humidities at ddraught top
    do jl = kidia , kfdia
      if ( llddraf(jl) ) then
        jk = idtop(jl)
        ik = min(jk+1,klev)
        if ( zmfdq(jl,jk) < 0.3D0*zmfdq(jl,ik) ) then
          if ( rmfsoltq == d_zero ) then
            zmfdq(jl,jk) = 0.3D0*zmfdq(jl,ik)
          else
            pmfd(jl,jk) = 0.3D0*pmfd(jl,ik)
          end if
        end if
      end if
    end do

    ! avoid negative humidities near cloud top because gradient of precip flux
    ! and detrainment / liquid water flux too large
    do jk = 2 , klev
      do jl = kidia , kfdia
        if ( ldcum(jl) .and. jk >= kctop(jl)-1 .and. jk < kcbot(jl) ) then
          zdz = ptsphy*egrav/(paph(jl,jk+1)-paph(jl,jk))
          zmfa = zmfuq(jl,jk+1)+zmfdq(jl,jk+1)-zmfuq(jl,jk)-zmfdq(jl,jk) + &
                 zmful(jl,jk+1)-zmful(jl,jk)+zdmfup(jl,jk)
          zmfa = (zmfa-plude(jl,jk))*zdz
          if ( pqen(jl,jk)+zmfa < d_zero ) then
            plude(jl,jk) = plude(jl,jk)+2.0D0*(pqen(jl,jk)+zmfa)/zdz
          end if
          if ( plude(jl,jk) < d_zero ) then
            plude(jl,jk) = d_zero
          end if
        end if
        if ( .not. ldcum(jl) ) then
          pmfude_rate(jl,jk) = d_zero
        end if
        if ( pmfd(jl,jk-1) == d_zero ) then
          pmfdde_rate(jl,jk) = d_zero
        end if
      end do
    end do

    !*upg change to operations
    if ( llconscheck ) then
      allocate(ztent(klon,klev))
      allocate(ztenq(klon,klev))
      do jk = 2 , klev
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            ztent(jl,jk) = ptent(jl,jk)
            ztenq(jl,jk) = ptenq(jl,jk)
            ztenu(jl,jk) = ptenu(jl,jk)
            ztenv(jl,jk) = ptenv(jl,jk)
          end if
        end do
      end do
      if ( lmftrac .and. ktrac > 0 ) then
        allocate(ztenc(klon,klev,ktrac))
        allocate(zsumc(klon,4+ktrac))
        do jn = 1 , ktrac
          do jk = 2 , klev
            do jl = kidia , kfdia
              if ( ldcum(jl) ) then
                ztenc(jl,jk,jn) = ptenc(jl,jk,jn)
              end if
            end do
          end do
        end do
      else
        allocate(zsumc(klon,4))
      end if
    end if
    !*upg change to operations

    !----------------------------------------------------------------------
    !*    8.0          update tendencies for t and q in subroutine cudtdq
    !                  --------------------------------------------------

    if ( rmfsoltq > d_zero ) then
      ! derive draught properties for implicit
      do jk = klev , 2 , -1
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            if ( jk > kcbot(jl) ) then
              zmfa = d_one/max(1.D-15,pmfu(jl,jk))
              pqu(jl,jk) = zqenh(jl,jk)+zmfuq(jl,jk)*zmfa
              ptu(jl,jk) = ztenh(jl,jk)+zmfus(jl,jk)*zmfa*rcpd
              zmfus(jl,jk) = pmfu(jl,jk)*(cpd*ptu(jl,jk)+pgeoh(jl,jk))
              zmfuq(jl,jk) = pmfu(jl,jk)*pqu(jl,jk)
              if ( llddraf(jl) ) then
                zmfa = d_one/min(-1.D-15,pmfd(jl,jk))
                zqd(jl,jk) = zqenh(jl,jk)+zmfdq(jl,jk)*zmfa
                ztd(jl,jk) = ztenh(jl,jk)+zmfds(jl,jk)*zmfa*rcpd
                zmfdq(jl,jk) = pmfd(jl,jk)*zqd(jl,jk)
                zmfds(jl,jk) = pmfd(jl,jk)*(cpd*ztd(jl,jk)+pgeoh(jl,jk))
              end if
            else if ( jk <= kcbot(jl) .and. jk >= kctop(jl) ) then
              zmfus(jl,jk) = pmfu(jl,jk)*(cpd*ptu(jl,jk)+pgeoh(jl,jk))
              zmfuq(jl,jk) = pmfu(jl,jk)*pqu(jl,jk)
              zmfds(jl,jk) = pmfd(jl,jk)*(cpd*ztd(jl,jk)+pgeoh(jl,jk))
              zmfdq(jl,jk) = pmfd(jl,jk)*zqd(jl,jk)
            end if
          end if
        end do
      end do
    end if

    call cudtdqn(kidia,kfdia,klon,klev,itopm2,kctop,idtop,ldcum,llddraf, &
                 ptsphy,paph,pgeoh,pgeo,pten,ztenh,pqen,zqenh,pqsen,     &
                 zlglac,plude,pmfu,pmfd,zmfus,zmfds,zmfuq,zmfdq,zmful,   &
                 zdmfup,zdpmel,ptent,ptenq,penth)

    !----------------------------------------------------------------------
    !*    9.0          compute momentum in updraught and downdraught
    !                  ---------------------------------------------

    if ( lmfdudv ) then

      do jk = klev-1 , 2 , -1
        ik = jk+1
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            if ( jk == kcbot(jl) .and. ktype(jl) < 3 ) then
              ikb = idpl(jl)
              zuu(jl,jk) = puen(jl,ikb-1)
              zvu(jl,jk) = pven(jl,ikb-1)
            else if ( jk == kcbot(jl) .and. ktype(jl) == 3 ) then
              zuu(jl,jk) = puen(jl,jk-1)
              zvu(jl,jk) = pven(jl,jk-1)
            end if
            if ( jk < kcbot(jl) .and. jk >= kctop(jl) ) then
              zfac = d_zero
              if ( ktype(jl) == 1 .or. ktype(jl) == 3 ) zfac = 2.0D0
              if ( ktype(jl) == 1 .and. jk <= kctop(jl)+2 ) zfac = 3.0D0
              zerate = pmfu(jl,jk)-pmfu(jl,ik)+(d_one+zfac)*pmfude_rate(jl,jk)
              zderate = (d_one+zfac)*pmfude_rate(jl,jk)
              zmfa = d_one/max(rmfcmin,pmfu(jl,jk))
              zuu(jl,jk) = (zuu(jl,ik)*pmfu(jl,ik) + &
                zerate*puen(jl,jk)-zderate*zuu(jl,ik))*zmfa
              zvu(jl,jk) = (zvu(jl,ik)*pmfu(jl,ik) + &
                zerate*pven(jl,jk)-zderate*zvu(jl,ik))*zmfa
            end if
          end if
        end do
      end do
      do jk = 3 , klev
        ik = jk-1
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            if ( jk == idtop(jl) ) then
              zud(jl,jk) = 0.5D0*(zuu(jl,jk)+puen(jl,ik))
              zvd(jl,jk) = 0.5D0*(zvu(jl,jk)+pven(jl,ik))
            else if ( jk > idtop(jl) ) then
              zerate = -pmfd(jl,jk)+pmfd(jl,ik)+pmfdde_rate(jl,jk)
              zmfa = d_one/min(-rmfcmin,pmfd(jl,jk))
              zud(jl,jk) = (zud(jl,ik)*pmfd(jl,ik) - &
                zerate*puen(jl,ik)+pmfdde_rate(jl,jk)*zud(jl,ik))*zmfa
              zvd(jl,jk) = (zvd(jl,ik)*pmfd(jl,ik) - &
                zerate*pven(jl,ik)+pmfdde_rate(jl,jk)*zvd(jl,ik))*zmfa
            end if
          end if
        end do
      end do

      !*    9.1          update tendencies for u and v in subroutine cududv
      !                  --------------------------------------------------
      ! for explicit/semi-implicit rescale massfluxes for stability in momentum
      !------------------------------------------------------------------------

      zmfs(:) = d_one
      ! if ( rmfsoluv <= 0.5D0 ) then
      if ( rmfsoluv <= d_one ) then
        do jk = 2 , klev
          do jl = kidia , kfdia
            if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
              zmfmax = (paph(jl,jk)-paph(jl,jk-1))*zcons
              if ( pmfu(jl,jk) > zmfmax .and. jk >= kctop(jl) ) then
                zmfs(jl) = min(zmfs(jl),zmfmax/pmfu(jl,jk))
              end if
            end if
          end do
        end do
      end if
      do jk = 1 , klev
        do jl = kidia , kfdia
          zmfuus(jl,jk) = pmfu(jl,jk)
          zmfdus(jl,jk) = pmfd(jl,jk)
          if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
            zmfuus(jl,jk) = pmfu(jl,jk)*zmfs(jl)
            zmfdus(jl,jk) = pmfd(jl,jk)*zmfs(jl)
          end if
        end do
      end do

      ! recompute draught properties below for implicit
      ! based on linear flux profiles

      if ( rmfsoluv > d_zero ) then
        do jl = kidia , kfdia
          if ( ldcum(jl) ) then
            jk = kcbot(jl)
            ik = jk-1
            zmfuub(jl) = zmfuus(jl,jk)*(zuu(jl,jk)-puen(jl,ik))
            zmfuvb(jl) = zmfuus(jl,jk)*(zvu(jl,jk)-pven(jl,ik))
          end if
        end do
        do jk = 2 , klev
          ik = jk-1
          do jl = kidia , kfdia
            if ( ldcum(jl) .and. jk > kcbot(jl) ) then
              ikb = kcbot(jl)
              zdz = ((paph(jl,klev+1)-paph(jl,jk)) / &
                     (paph(jl,klev+1)-paph(jl,ikb)))
              if ( ktype(jl) == 3 ) then
                zdz = zdz*zdz
              end if
              zmfa = d_one/max(rmfcmin,zmfuus(jl,jk))
              zuu(jl,jk) = puen(jl,ik)+zmfuub(jl)*zdz*zmfa
              zvu(jl,jk) = pven(jl,ik)+zmfuvb(jl)*zdz*zmfa

              zmfdus(jl,jk) = zmfdus(jl,ikb)*zdz
              zud(jl,jk) = puen(jl,ik)+zud(jl,ikb)-puen(jl,ikb-1)
              zvd(jl,jk) = pven(jl,ik)+zvd(jl,ikb)-pven(jl,ikb-1)
            end if
            ! add uv perturb to correct wind bias
            if ( ldcum(jl) .and. jk >= kctop(jl) ) then
              zuu(jl,jk) = zuu(jl,jk)-ruvper*sign(d_one,zuu(jl,jk))
              zvu(jl,jk) = zvu(jl,jk)-ruvper*sign(d_one,zvu(jl,jk))
            end if
          end do
        end do
      end if

      !-------------------------------------------------------------------
      ! end
      ! intermediate solution for stability in eps:
      ! for original code replace line
      !  &, puen,     pven,     zmfuus,   zmfdus &
      !by
      !  &, puen,     pven,     pmfu,     pmfd

      call cududv(kidia,kfdia,klon,klev,itopm2,ktype,kcbot,kctop, &
                  ldcum,ptsphy,paph,puen,pven,zmfuus,zmfdus,zuu,  &
                  zud,zvu,zvd,ptenu,ptenv)

      if ( lmfuvdis ) then
        ! add ke dissipation
        do jl = kidia , kfdia
          zsum12(jl) = d_zero
          zsum22(jl) = d_zero
        end do
        do jk = 1 , klev
          do jl = kidia , kfdia
            zuv2(jl,jk) = d_zero
            if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
              zdz = (paph(jl,jk+1)-paph(jl,jk))
              zduten = ptenu(jl,jk)-ztenu(jl,jk)
              zdvten = ptenv(jl,jk)-ztenv(jl,jk)
              zuv2(jl,jk) = sqrt(zduten**2+zdvten**2)
              zsum22(jl) = zsum22(jl)+zuv2(jl,jk)*zdz
              zsum12(jl) = zsum12(jl) - &
                (puen(jl,jk)*zduten+pven(jl,jk)*zdvten)*zdz
            end if
          end do
        end do
        do jk = 1 , klev
          do jl = kidia , kfdia
            if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
              zdz = (paph(jl,jk+1)-paph(jl,jk))
              ztdis = rcpd*zsum12(jl)*zuv2(jl,jk)/max(1.D-15,zsum22(jl))
              ptent(jl,jk) = ptent(jl,jk)+ztdis
            end if
          end do
        end do
      end if
    end if

    !----------------------------------------------------------------------
    !*   10.           in case that either deep or shallow is switched off
    !                  need to set some variables a posteriori to zero
    !                  ---------------------------------------------------

    if ( .not. lmfscv .or. .not. lmfpen ) then
      do jk = 2 , klev
        do jl = kidia , kfdia
          if ( llo2(jl) .and. jk >= kctop(jl)-1 ) then
            ptu(jl,jk) = pten(jl,jk)
            pqu(jl,jk) = pqen(jl,jk)
            plu(jl,jk) = d_zero
            penth(jl,jk) = d_zero
            pmfude_rate(jl,jk) = d_zero
            pmfdde_rate(jl,jk) = d_zero
          end if
        end do
      end do
      do jl = kidia , kfdia
        if ( llo2(jl) ) then
          kctop(jl) = klev-1
          kcbot(jl) = klev-1
        end if
      end do
    end if

    !----------------------------------------------------------------------
    !*   11.0          chemical tracer transport
    !                  -------------------------

    if ( lmftrac .and. ktrac > 0 ) then

      ! transport switched off for mid-level convection
      do jl = kidia , kfdia
        ! if ( ldcum(jl) .and. ktype(jl) /= 3 ) then
         if ( ldcum(jl) .and. ktype(jl) /= 3 .and. &
              kcbot(jl)-kctop(jl) >= 1 ) then
           lldcum(jl) = .true.
           llddraf3(jl) = llddraf(jl)
         else
           lldcum(jl) = .false.
           llddraf3(jl) = .false.
         end if
      end do

      ! check and correct mass fluxes for cfl criterium

      zmfs(:) = d_one
      if ( rmfsolct <= 3.0D0 ) then
        do jk = 2 , klev
          do jl = kidia , kfdia
            if ( lldcum(jl) .and. jk >= kctop(jl) ) then
              zmfmax = (paph(jl,jk)-paph(jl,jk-1))*0.8D0*zcons
              if ( pmfu(jl,jk) > zmfmax ) then
                zmfs(jl) = min(zmfs(jl),zmfmax/pmfu(jl,jk))
              end if
            end if
          end do
        end do
      end if
      do jk = 1 , klev
        do jl = kidia , kfdia
          if ( lldcum(jl) .and. jk >= kctop(jl)-1 ) then
            zmfuus(jl,jk) = pmfu(jl,jk)*zmfs(jl)
            zmfudr(jl,jk) = pmfude_rate(jl,jk)*zmfs(jl)
          else
            zmfuus(jl,jk) = d_zero
            zmfudr(jl,jk) = d_zero
          end if
          if ( llddraf3(jl) .and. jk >= idtop(jl)-1 ) then
            zmfdus(jl,jk) = pmfd(jl,jk)*zmfs(jl)
            zmfddr(jl,jk) = pmfdde_rate(jl,jk)*zmfs(jl)
          else
            zmfdus(jl,jk) = d_zero
            zmfddr(jl,jk) = d_zero
          end if
        end do
      end do

      if ( lmfsmooth ) then
        ! smoothing of mass fluxes (gradients) at top and bottom of draughts
        do jk = 2 , klev-1
          do jl = kidia , kfdia
            if ( llddraf3(jl) .and. zmfdus(jl,jk) < d_zero .and. &
                 zmfdus(jl,jk+1) == d_zero ) then
              zerate = min(d_zero,zmfdus(jl,jk)-0.5*zmfdus(jl,jk-1))
              zmfdus(jl,jk) = zmfdus(jl,jk)-zerate
              zmfddr(jl,jk) = zmfddr(jl,jk)-zerate
              zmfddr(jl,jk+1) = -zmfdus(jl,jk)
            end if
            if ( lldcum(jl) .and. jk == kctop(jl) ) then
              zerate = max(d_zero,zmfuus(jl,jk)-0.5D0*zmfuus(jl,jk+1))
              zmfuus(jl,jk) = zmfuus(jl,jk)-zerate
              zmfudr(jl,jk) = zmfudr(jl,jk)+zerate
              zmfudr(jl,jk-1) = zmfuus(jl,jk)
            end if
          end do
        end do
        do jk = klev-1 , 2 , -1
          do jl = kidia , kfdia
            if ( lldcum(jl) ) then
              if ( zmfudr(jl,jk) == d_zero .and. &
                   zmfudr(jl,jk-1) > d_zero ) then
                zmfudr(jl,jk) = 0.5D0*zmfudr(jl,jk-1)
              end if
            end if
          end do
        end do
      end if

      call cuctracer(kidia,kfdia,klon,ktdia,klev,ktrac,ktype,kctop,kcbot, &
                     idpl,idtop,lldcum,llddraf3,ptsphy,paph,zmfuus,zmfdus,&
                     zmfudr,zmfddr,pcen,ptenc)
    end if

    !----------------------------------------------------------------------
    !*   12.           put detrainment rates from mflx units in units mflx/m
    !                  for era40
    !                  ---------------------------------------------------

    do jk = 2 , klev
      do jl = kidia , kfdia
        if ( ldcum(jl) ) then
          zro = egrav/(pgeoh(jl,jk)-pgeoh(jl,jk+1))  ! 1/dz
          pmfude_rate(jl,jk) = pmfude_rate(jl,jk)*zro
          pmfdde_rate(jl,jk) = pmfdde_rate(jl,jk)*zro
          if ( jk<kctop(jl) ) then
            plu(jl,jk) = d_zero
            ptu(jl,jk) = pten(jl,jk)
            pqu(jl,jk) = pqen(jl,jk)
          end if
        end if
      end do
    end do

    !----------------------------------------------------------------------
    !*upg change to operations

    if ( llconscheck ) then

      !*   13.0          conservation check and correction
      !                  ---------------------------------

      do jl = kidia , kfdia
        zsumc(jl,:) = d_zero
      end do
      do jk = klev , 2 , -1
        do jl = kidia , kfdia
          if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
            zdz = (paph(jl,jk+1)-paph(jl,jk))/egrav
            zsumc(jl,1) = zsumc(jl,1) + &
              (ptenq(jl,jk)-ztenq(jl,jk))*zdz+plude(jl,jk)
            zalv = foelhmcu(pten(jl,jk))
            zsumc(jl,2) = zsumc(jl,2) + &
              cpd*(ptent(jl,jk)-ztent(jl,jk))*zdz-zalv*plude(jl,jk)
            zsumc(jl,3) = zsumc(jl,3)+(ptenu(jl,jk)-ztenu(jl,jk))*zdz
            zsumc(jl,4) = zsumc(jl,4)+(ptenv(jl,jk)-ztenv(jl,jk))*zdz
          end if
        end do
      end do
      if ( lmftrac .and. ktrac > 0 ) then
        do jn = 1 , ktrac
          do jk = klev , 2 , -1
            do jl = kidia , kfdia
              if ( ldcum(jl) .and. jk >= kctop(jl)-1 ) then
                zdz = (paph(jl,jk+1)-paph(jl,jk))/egrav
                zsumc(jl,4+jn) = zsumc(jl,4+jn) + &
                  (ptenc(jl,jk,jn)-ztenc(jl,jk,jn))*zdz
              end if
            end do
          end do
        end do
      end if

      do jl = kidia , kfdia
        if ( ldcum(jl) ) then
          zalv = foelhmcu(pten(jl,klev))
          zsfl(jl) = pmflxr(jl,klev+1)+pmflxs(jl,klev+1)
          write(61,'(i4,a9,2f15.8,i4,a9,f15.8,a10,2f15.8)')              &
            jl,' cons q: ',-zsumc(jl,1)*zalv,zsfl(jl)*zalv,ktype(jl),    &
            ' cons h: ',zsumc(jl,2),' cons uv: ',zsumc(jl,3),zsumc(jl,4)
          if ( lmftrac .and. ktrac > 0 ) then
            write(61,*)' conserv error tracers 1-',ktrac,' :'
            do jn = 1 , ktrac
              write(61,'(i4,e12.4)') jn,zsumc(jl,4+jn)
            end do
          end if

          ikb = kctop(jl)
          zdz = (paph(jl,klev+1)-paph(jl,ikb-1))/egrav
          zsumc(jl,1) = (zsumc(jl,1)+zsfl(jl))/zdz
          zsumc(jl,2) = (zsumc(jl,2)-zalv*zsfl(jl))/(zdz*cpd)
        end if
      end do

      ! do jk = klev , 2 , -1
      !   do jl = kidia , kfdia
      !     ikb = kctop(jl)
      !     if ( ldcum(jl) .and. jk >= ikb-1 ) then
      !       ptenq(jl,jk) = ptenq(jl,jk)-zsumc(jl,1)
      !       ptent(jl,jk) = ptent(jl,jk)-zsumc(jl,2)
      !       ptenu(jl,jk) = ptenu(jl,jk)-zsumc(jl,3)
      !       ptenv(jl,jk) = ptenv(jl,jk)-zsumc(jl,4)
      !     end if
      !   end do
      ! end do

      deallocate(zsumc)
      if ( lmftrac .and. ktrac > 0 ) then
        deallocate(ztenc)
      end if
      deallocate(ztenq)
      deallocate(ztent)

    end if

    !----------------------------------------------------------------------
    !*    14.0         compute convective tendencies for liquid and solid
    !                  cloud condensate, change precip units in m/s
    !                  --------------------------------------------------

    do jk = 1 , klev
      do jl = kidia , kfdia
        ptenl(jl,jk) = plude(jl,jk)*egrav/(paph(jl,jk+1)-paph(jl,jk))
        pteni(jl,jk) = (d_one-foealfa(pten(jl,jk)))*ptenl(jl,jk)
        ptenl(jl,jk) = ptenl(jl,jk)-pteni(jl,jk)
        pmflxr(jl,jk) = pmflxr(jl,jk)*1.D-3
        pmflxs(jl,jk) = pmflxs(jl,jk)*1.D-3
      end do
    end do
    do jl = kidia , kfdia
      pmflxr(jl,klev+1) = pmflxr(jl,klev+1)*1.D-3
      pmflxs(jl,klev+1) = pmflxs(jl,klev+1)*1.D-3
    end do
    !----------------------------------------------------------------------
    !*upg change to operations

  end subroutine cumastrn
!
!-----------------------------------------------------------------------------
!
  real(rk8) function minj(x,y)
    implicit none
    real(rk8) , intent(in) :: x , y
    minj = y - d_half*(dabs(x-y)-(x-y))
  end function minj
  real(rk8) function maxj(x,y)
    implicit none
    real(rk8) , intent(in) :: x , y
    maxj = y + d_half*(dabs(x-y)+(x-y))
  end function maxj
  real(rk8) function foedelta(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foedelta = max(d_zero,sign(d_one,ptare-tzero))
  end function foedelta
  real(rk8) function foeew(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foeew = c2es*exp((c3les*foedelta(ptare) + &
            c3ies*(d_one-foedelta(ptare)))*(ptare-tzero) / &
            (ptare-(c4les*foedelta(ptare)+c4ies*(d_one-foedelta(ptare)))))
  end function foeew
  real(rk8) function foede(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foede = (foedelta(ptare)*c5alvcp+(d_one-foedelta(ptare))*c5alscp) / &
       (ptare-(c4les*foedelta(ptare)+c4ies*(d_one-foedelta(ptare))))**2
  end function foede
  real(rk8) function foedesu(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foedesu = (foedelta(ptare)*c5les+(d_one-foedelta(ptare))*c5ies) / &
         (ptare-(c4les*foedelta(ptare)+c4ies*(d_one-foedelta(ptare))))**2
  end function foedesu
  real(rk8) function foelh(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foelh = foedelta(ptare)*wlhv + (d_one-foedelta(ptare))*wlhs
  end function foelh
  real(rk8) function foeldcp(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foeldcp = foedelta(ptare)*wlhvocp + (d_one-foedelta(ptare))*wlhsocp
  end function foeldcp
  real(rk8) function foealfa(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foealfa = min(d_one,((max(rtice,min(rtwat,ptare))-rtice) * &
                  rtwat_rtice_r)**2)
  end function foealfa
  real(rk8) function foeewm(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foeewm = c2es*(foealfa(ptare)*exp(c3les*(ptare-tzero)/(ptare-c4les))+ &
          (d_one-foealfa(ptare))*exp(c3ies*(ptare-tzero)/(ptare-c4ies)))
  end function foeewm
  real(rk8) function foedem(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foedem = foealfa(ptare)*c5alvcp*(d_one/(ptare-c4les)**2) + &
            (d_one-foealfa(ptare))*c5alscp*(d_one/(ptare-c4ies)**2)
  end function foedem
  real(rk8) function foeldcpm(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foeldcpm = foealfa(ptare)*wlhvocp+(d_one-foealfa(ptare))*wlhsocp
  end function foeldcpm
  real(rk8) function foelhm(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foelhm = foealfa(ptare)*wlhv+(d_one-foealfa(ptare))*wlhs
  end function foelhm
  real(rk8) function foetb(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foetb = foealfa(ptare)*c3les*(tzero-c4les)*(d_one/(ptare-c4les)**2)+ &
      (d_one-foealfa(ptare))*c3ies*(tzero-c4ies)*(d_one/(ptare-c4ies)**2)
  end function foetb
  real(rk8) function foealfcu(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foealfcu = min(d_one, &
           ((max(rtice,min(rtwat,ptare))-rtice)*rtwat_rtice_r)**2)
  end function foealfcu
  real(rk8) function foeewmcu(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foeewmcu = c2es*(foealfcu(ptare)*exp(c3les*(ptare-tzero)/(ptare-c4les))+ &
            (d_one-foealfcu(ptare))*exp(c3ies*(ptare-tzero)/(ptare-c4ies)))
  end function foeewmcu
  real(rk8) function foedemcu(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foedemcu = foealfcu(ptare)*c5alvcp*(d_one/(ptare-c4les)**2) + &
           (d_one-foealfcu(ptare))*c5alscp*(d_one/(ptare-c4ies)**2)
  end function foedemcu
  real(rk8) function foeldcpmcu(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foeldcpmcu = foealfcu(ptare)*wlhvocp+(d_one-foealfcu(ptare))*wlhsocp
  end function foeldcpmcu
  real(rk8) function foelhmcu(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foelhmcu = foealfcu(ptare)*wlhv+(d_one-foealfcu(ptare))*wlhs
  end function foelhmcu
  real(rk8) function foeewmo(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foeewmo = c2es*exp(c3les*(ptare-tzero)/(ptare-c4les))
  end function foeewmo
  real(rk8) function foeeliq(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foeeliq = c2es*exp(c3les*(ptare-tzero)/(ptare-c4les))
  end function foeeliq
  real(rk8) function foeeice(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foeeice = c2es*exp(c3ies*(ptare-tzero)/(ptare-c4ies))
  end function foeeice
  real(rk8) function foeles_v(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foeles_v = c3les*(ptare-tzero)/(ptare-c4les)
  end function foeles_v
  real(rk8) function foeies_v(ptare)
    implicit none
    real(rk8) , intent(in) :: ptare
    foeies_v = c3ies*(ptare-tzero)/(ptare-c4ies)
  end function foeies_v
  real(rk8) function foeewm_v(ptare,exp1,exp2)
    implicit none
    real(rk8) , intent(in) :: ptare , exp1 , exp2
    foeewm_v = c2es*(foealfa(ptare)*exp1+(d_one-foealfa(ptare))*exp2)
  end function foeewm_v
  real(rk8) function foeewmcu_v(ptare,exp1,exp2)
    implicit none
    real(rk8) , intent(in) :: ptare , exp1 , exp2
    foeewmcu_v = c2es*(foealfcu(ptare)*exp1+(d_one-foealfcu(ptare))*exp2)
  end function foeewmcu_v
  subroutine vdiv(z,x,y,n)
    implicit none
    real(rk8) , dimension(:) , intent(out) :: z
    real(rk8) , dimension(:) , intent(in) :: x , y
    integer(ik4) , intent(in) :: n
    integer(ik4) :: i
    do i = 1 , n
      z(i) = x(i)/y(i)
    end do
  end subroutine vdiv
  subroutine vexp(y,x,n)
    implicit none
    real(rk8) , dimension(:) , intent(out) :: y
    real(rk8) , dimension(:) , intent(in) :: x
    integer(ik4) , intent(in) :: n
    integer(ik4) :: i
    do i = 1 , n
      y(i) = dexp(x(i))
    end do
  end subroutine vexp
  subroutine vrec(y,x,n)
    implicit none
    real(rk8) , dimension(:) , intent(out) :: y
    real(rk8) , dimension(:) , intent(in) :: x
    integer(ik4) , intent(in) :: n
    integer(ik4) :: i
    do i = 1 , n
      y(i) = d_one/x(i)
    end do
  end subroutine vrec
end module mod_cu_tiedtke_38r2
