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

module mod_cloud_thomp

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams

  implicit none

  private

  public :: thomp_cldfrac

  contains
  !
  ! Cloud fraction scheme by G. Thompson (NCAR-RAL), not intended for
  ! combining with any cumulus or shallow cumulus parameterization
  ! scheme cloud fractions.  This is intended as a stand-alone for
  ! cloud fraction and is relatively good at getting widespread stratus
  ! and stratoCu without caring whether any deep/shallow Cu param schemes
  ! is making sub-grid-spacing clouds/precip.  Under the hood, this
  ! scheme follows Mocko and Cotton (1995) in applicaiton of the
  ! Sundqvist et al (1989) scheme but using a grid-scale dependent
  ! RH threshold, one each for land v. ocean points.
  !
  subroutine thomp_cldfrac(p,t,qv,qc,qs,qi,iland,gridkm,cldfra)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: p, t
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: qv, qc, qi, qs
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: iland
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: cldfra
    real(rkx), intent(in) :: gridkm
    real(rkx) :: rh_00l, rh_00o
    real(rkx), dimension(jci1:jci2,ici1:ici2,1:kz):: qvsat
    integer(ik4) :: i, j, k
    real(rkx) :: rh_00, rhi_max
    real(rkx) :: tk, tc, qvsi, qvsw, rhum

    ! First cut scale-aware. Higher resolution should require closer to
    ! saturated grid box for higher cloud fraction.  Simple functions
    ! chosen based on mocko and cotton (1995) starting point and desire
    ! to get near 100% rh as grid spacing moves toward 1.0km, but higher
    ! rh over ocean required as compared to over land.

    rh_00l = 0.781_rkx + sqrt(d_one/(50.0_rkx+gridkm*gridkm*gridkm*0.5_rkx))
    rh_00o = 0.831_rkx + sqrt(d_one/(70.0_rkx+gridkm*gridkm*gridkm*0.5_rkx))

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      cldfra(j,i,k) = d_zero
      rhi_max = d_zero
      if ( qc(j,i,k) > 1.e-4_rkx .or.  &
           qi(j,i,k) >= 1.e-6_rkx .or. &
           qs(j,i,k) > 1.e-4_rkx ) then
        cldfra(j,i,k) = d_one
        qvsat(j,i,k) = qv(j,i,k)
      else
        tk = t(j,i,k)
        tc = tk - tzero
        qvsw = rslf(p(j,i,k), tk)
        qvsi = rsif(p(j,i,k), tk)
        if ( tc >= -12.0_rkx ) then
          qvsat(j,i,k) = qvsw
        else if ( tc < -20.0_rkx )  then
          qvsat(j,i,k) = qvsi
        else
          qvsat(j,i,k) = qvsw - &
            (qvsw-qvsi)*(-12.0_rkx-tc)/(-12.0_rkx+20.0_rkx)
        end if
        rhum = max(0.01_rkx,min(0.999_rkx,qv(j,i,k)/qvsat(j,i,k)))
        if ( iland(j,i) == 0 ) then
          rh_00 = rh_00o
        else
          rh_00 = rh_00l
        end if
        if ( tc >= -12.0_rkx ) then
          cldfra(j,i,k) = max(d_zero,d_one - &
                      sqrt((d_one-rhum)/(d_one-rh_00)))
        else if ( tc < -12.0_rkx .and. tc > -70.0_rkx .and. &
                  rhum > rh_00o ) then
          rhi_max = max(rhum+1.e-6_rkx, qvsw/qvsi)
          cldfra(j,i,k) = max(d_zero,d_one - &
                      sqrt((rhi_max-rhum)/(rhi_max-rh_00)))
        end if
        cldfra(j,i,k) = max(d_zero, min(cldfra(j,i,k), d_one))
      end if
    end do

    contains
     !
     ! This function calculates the liquid saturation vapor mixing ratio as
     ! a function of temperature and pressure
     !
     pure real(rkx) function rslf(p,t)
!$acc routine seq
       implicit none
       real(rkx), intent(in) :: p, t
       rslf = pfwsatl(t,p)
     end function rslf
     !
     ! This function calculates the ice saturation vapor mixing ratio as a
     ! function of temperature and pressure
     !
     pure real(rkx) function rsif(p,t)
!$acc routine seq
       implicit none
       real(rkx), intent(in) :: p, t
       real(rkx) :: esi
       esi = min(esatice(t),0.15_rkx*p)
       rsif = ep2 * esi/(p-esi)
     end function rsif

  end subroutine thomp_cldfrac

#include <saturation.inc>

end module mod_cloud_thomp
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
