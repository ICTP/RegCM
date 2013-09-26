      MODULE mod_che_common_isorropia

!=======================================================================
! *** ISORROPIA CODE II
! *** INCLUDE FILE 'ISRPIA.INC'
! *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
!     AND VARIABLES. 
!
! *** COPYRIGHT 1996-2008, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS
!
!=======================================================================

!      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      ! INTEGER parameters
      INTEGER, PARAMETER :: NCOMP  =   8
      INTEGER, PARAMETER :: NIONS  =  10
      INTEGER, PARAMETER :: NGASAQ =   3
      INTEGER, PARAMETER :: NSLDS  =  19
      INTEGER, PARAMETER :: NPAIR  =  23
      INTEGER, PARAMETER :: NZSR   = 100
      INTEGER, PARAMETER :: NERRMX =  25

      !=================================================================
      ! Input variables
      !=================================================================

      INTEGER            :: IPROB, METSTBL
      INTEGER            :: NADJ
      REAL*8             :: W(NCOMP), WAER(NCOMP), TEMP, RH
!      COMMON /INPT/ W, WAER, TEMP, RH, IPROB, METSTBL, NADJ

      !=================================================================
      ! Water activities of pure salt solutions
      !=================================================================

      REAL*8     :: AWAS(NZSR), AWSS(NZSR), AWAC(NZSR), AWSC(NZSR)
      REAL*8     :: AWAN(NZSR), AWSN(NZSR), AWSB(NZSR), AWAB(NZSR)
      REAL*8     :: AWSA(NZSR), AWLC(NZSR), AWCS(NZSR), AWCN(NZSR)
      REAL*8     :: AWCC(NZSR), AWPS(NZSR), AWPB(NZSR), AWPN(NZSR)
      REAL*8     :: AWPC(NZSR), AWMS(NZSR), AWMN(NZSR), AWMC(NZSR)
!      COMMON /ZSR / AWAS,       AWSS,       AWAC,       AWSC,
!     &              AWAN,       AWSN,       AWSB,       AWAB,
!     &              AWSA,       AWLC,       AWCS,       AWCN,
!     &              AWCC,       AWPS,       AWPB,       AWPN,
!     &              AWPC,       AWMS,       AWMN,       AWMC

      !=================================================================
      ! Deliquescence relative humidities
      !=================================================================

      INTEGER WFTYP
      REAL*8     :: DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3 
      REAL*8     :: DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL
      REAL*8     :: DRCASO4,  DRCANO32, DRCACL2,  DRK2SO4,  DRKHSO4
      REAL*8     :: DRKNO3,   DRKCL,    DRMGSO4,  DRMGNO32, DRMGCL2
!      COMMON /DRH / DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3, 
!     &              DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL,
!     &              DRCASO4,  DRCANO32, DRCACL2,  DRK2SO4,  DRKHSO4,
!     &              DRKNO3,   DRKCL,    DRMGSO4,  DRMGNO32, DRMGCL2

      REAL*8     ::  DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2
      REAL*8     ::  DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2
      REAL*8     ::  DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3
      REAL*8     ::  DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8
      REAL*8     ::  DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13
!      COMMON /MDRH/ DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2,
!     &              DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2,
!     &              DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3,
!     &              DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8,
!     &              DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13,
!     &              WFTYP

      REAL*8     ::  DRMO1,    DRMO2,    DRMO3,    DRML1,    DRML2
      REAL*8     ::  DRML3,    DRMM1,    DRMM2,    DRMP1,    DRMP2
      REAL*8     ::  DRMP3,    DRMP4,    DRMP5,    DRMV1
!      COMMON /MDRH2/ DRMO1,    DRMO2,    DRMO3,    DRML1,    DRML2,
!     &               DRML3,    DRMM1,    DRMM2,    DRMP1,    DRMP2,
!     &               DRMP3,    DRMP4,    DRMP5,    DRMV1

      !=================================================================
      ! Variables for liquid aerosol phase
      !=================================================================

      REAL*8 MOLAL(NIONS), MOLALR(NPAIR), M0(NPAIR)
      REAL*8 ZZ(NPAIR), Z(NIONS), EPSACT
      REAL*8 GAMOU(NPAIR), GAMIN(NPAIR), GASAQ(NGASAQ)
      REAL*8 GAMA(NPAIR)
      REAL*8 COH, CHNO3, CHCL, WATER
      INTEGER :: IACALC
      REAL IONIC
      LOGICAL CALAOU, CALAIN, FRST, DRYF
!      COMMON /IONS/ MOLAL,        MOLALR,        GAMA,        ZZ,
!     &              Z,            GAMOU,         GAMIN,       M0, 
!     &              GASAQ,        EPSACT,        COH,         CHNO3,
!     &              CHCL,         WATER,         IONIC,       IACALC,
!     &              FRST,         CALAIN,        CALAOU,      DRYF
      
      !=================================================================
      ! Variables for solid aerosol phase
      !=================================================================

      REAL*8     ::  CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4 
      REAL*8     ::  CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC, CCASO4
      REAL*8     ::  CCANO32, CCACL2,  CK2SO4,  CKHSO4,  CKNO3, CKCL
      REAL*8     ::  CMGSO4,  CMGNO32, CMGCL2
!      COMMON /SALT/ CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4, 
!     &              CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC, CCASO4,
!     &              CCANO32, CCACL2,  CK2SO4,  CKHSO4,  CKNO3, CKCL,
!     &              CMGSO4,  CMGNO32, CMGCL2

      !=================================================================
      ! Variables for gas phase
      !=================================================================

      REAL*8     ::  GNH3, GHNO3, GHCL 
!      COMMON /GAS / GNH3, GHNO3, GHCL 

      !=================================================================
      ! Equilibrium constants
      !=================================================================

      REAL*8     :: XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10
      REAL*8     :: XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41
      REAL*8     :: XK42,XK15,XK16,XK17,XK18,XK19,XK20,XK23
      REAL*8     :: XK24,XK25
!      COMMON /EQUK/ XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10,
!     &              XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41,
!     &              XK42, XK15, XK16, XK17, XK18, XK19, XK20, XK23,
!     &              XK24, XK25
!     &              , XK26, XK27

      !=================================================================
      ! Molecular Weights
      !=================================================================

      REAL*8     :: R, IMW(NIONS), WMW(NCOMP), SMW(NPAIR)
!      COMMON /OTHR/ R, IMW,        WMW,        SMW

      !=================================================================
      ! Solution/info variables
      !=================================================================

      CHARACTER(LEN=15) :: SCASE
      REAL*8            :: SULRATW, SULRAT,  SODRAT,
     &                     SO4RAT,  CRNARAT, CRRAT
!      COMMON /CASE/ SULRATW, SULRAT, SODRAT, SO4RAT, CRNARAT, CRRAT,
!     &              SCASE

      REAL*8            :: EPS
      INTEGER           :: MAXIT, NSWEEP, NDIV, ICLACT
!      COMMON /SOLN/ EPS, MAXIT, NSWEEP, NDIV, ICLACT

      !=================================================================
      ! Error system
      !=================================================================

      CHARACTER(LEN=40) :: ERRMSG(NERRMX)
      INTEGER   ERRSTK(NERRMX), NOFER   
      LOGICAL   STKOFL 
!      COMMON /EROR/ STKOFL, NOFER, ERRSTK, ERRMSG

      !=================================================================
      ! Generic Variables
      !=================================================================

      CHARACTER(LEN=15) :: VERSION
      REAL*8            :: GREAT, TINY, TINY2, ZERO, ONE
!      COMMON /CGEN/ GREAT, TINY, TINY2, ZERO, ONE, VERSION

!
! *** END OF INCLUDE FILE **********************************************
!

      END MODULE mod_che_common_isorropia


      MODULE mod_che_common_SOLUT

      REAL*8      ::  CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8
      REAL*8      ::  CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15
      REAL*8      ::  CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6
      REAL*8      ::  PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13
      REAL*8      ::  PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6
      REAL*8      ::  A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17

      END MODULE mod_che_common_SOLUT


      MODULE mod_che_common_CASEG

      REAL*8      ::  CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, LAMDA
      REAL*8      ::  PSI1, PSI2, PSI3, PSI4, PSI5, PSI6, PSI7
      REAL*8      ::  A1,   A2,   A3,   A4,   A5,   A6,   A7

      END MODULE mod_che_common_CASEG


      MODULE mod_che_common_CASEJ

      REAL*8      ::  CHI1, CHI2, CHI3, LAMDA, KAPA, PSI1, PSI2, PSI3
      REAL*8      ::  A1,   A2,   A3
 
      END MODULE mod_che_common_CASEJ


      MODULE mod_che_common_CASEK

      REAL*8      ::  CHI1,CHI2,CHI3,CHI4,LAMDA,KAPA,PSI1,PSI2,PSI3
      REAL*8      ::  A1,   A2,   A3,   A4

      END MODULE mod_che_common_CASEK


      MODULE mod_che_common_CASEO

      REAL*8      ::  CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8
      REAL*8      ::  CHI9, LAMDA, PSI1, PSI2, PSI3, PSI4, PSI5
      REAL*8      ::  PSI6, PSI7, PSI8, PSI9,  A1,  A2,  A3,  A4
      REAL*8      ::  A5, A6, A7, A8, A9

      END MODULE mod_che_common_CASEO
