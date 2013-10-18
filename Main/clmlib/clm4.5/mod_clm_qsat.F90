
module mod_clm_qsat

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: QSatMod
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: QSat
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: QSat
!
! !INTERFACE:
  subroutine QSat (T, p, es, esdT, qs, qsdT)
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.
! Reference:  Polynomial approximations from:
!             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
!             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
!
! !USES:
    use mod_realkinds
    use mod_constants
!
! !ARGUMENTS:
    implicit none
    real(rk8), intent(in)  :: T        ! temperature (K)
    real(rk8), intent(in)  :: p        ! surface atmospheric pressure (pa)
    real(rk8), intent(out) :: es       ! vapor pressure (pa)
    real(rk8), intent(out) :: esdT     ! d(es)/d(T)
    real(rk8), intent(out) :: qs       ! humidity (kg/kg)
    real(rk8), intent(out) :: qsdT     ! d(qs)/d(T)
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine CanopyFluxesMod CanopyFluxesMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!
! !LOCAL VARIABLES:
!EOP
!
    real(rk8) :: T_limit
    real(rk8) :: td,vp,vp1,vp2
!
! For water vapor (temperature range 0C-100C)
!
    real(rk8), parameter :: a0 =  6.11213476D0
    real(rk8), parameter :: a1 =  0.444007856D0
    real(rk8), parameter :: a2 =  0.143064234D-01
    real(rk8), parameter :: a3 =  0.264461437D-03
    real(rk8), parameter :: a4 =  0.305903558D-05
    real(rk8), parameter :: a5 =  0.196237241D-07
    real(rk8), parameter :: a6 =  0.892344772D-10
    real(rk8), parameter :: a7 = -0.373208410D-12
    real(rk8), parameter :: a8 =  0.209339997D-15
!
! For derivative:water vapor
!
    real(rk8), parameter :: b0 =  0.444017302D0
    real(rk8), parameter :: b1 =  0.286064092D-01
    real(rk8), parameter :: b2 =  0.794683137D-03
    real(rk8), parameter :: b3 =  0.121211669D-04
    real(rk8), parameter :: b4 =  0.103354611D-06
    real(rk8), parameter :: b5 =  0.404125005D-09
    real(rk8), parameter :: b6 = -0.788037859D-12
    real(rk8), parameter :: b7 = -0.114596802D-13
    real(rk8), parameter :: b8 =  0.381294516D-16
!
! For ice (temperature range -75C-0C)
!
    real(rk8), parameter :: c0 =  6.11123516D0
    real(rk8), parameter :: c1 =  0.503109514D0
    real(rk8), parameter :: c2 =  0.188369801D-01
    real(rk8), parameter :: c3 =  0.420547422D-03
    real(rk8), parameter :: c4 =  0.614396778D-05
    real(rk8), parameter :: c5 =  0.602780717D-07
    real(rk8), parameter :: c6 =  0.387940929D-09
    real(rk8), parameter :: c7 =  0.149436277D-11
    real(rk8), parameter :: c8 =  0.262655803D-14
!
! For derivative:ice
!
    real(rk8), parameter :: d0 =  0.503277922D0
    real(rk8), parameter :: d1 =  0.377289173D-01
    real(rk8), parameter :: d2 =  0.126801703D-02
    real(rk8), parameter :: d3 =  0.249468427D-04
    real(rk8), parameter :: d4 =  0.313703411D-06
    real(rk8), parameter :: d5 =  0.257180651D-08
    real(rk8), parameter :: d6 =  0.133268878D-10
    real(rk8), parameter :: d7 =  0.394116744D-13
    real(rk8), parameter :: d8 =  0.498070196D-16
!-----------------------------------------------------------------------

    T_limit = T - tzero
    if (T_limit > 100.0D0) T_limit=100.0D0
    if (T_limit < -75.0D0) T_limit=-75.0D0

    td       = T_limit
    if (td >= 0.0D0) then
       es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
            + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
       esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
            + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
       es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
            + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
       esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
            + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    endif

    es    = es    * 100.D0            ! pa
    esdT  = esdT  * 100.D0            ! pa/K

    vp    = 1.0D0   / (p - 0.378D0*es)
    vp1   = 0.622D0 * vp
    vp2   = vp1   * vp

    qs    = es    * vp1             ! kg/kg
    qsdT  = esdT  * vp2 * p         ! 1 / K

  end subroutine QSat

end module mod_clm_qsat
