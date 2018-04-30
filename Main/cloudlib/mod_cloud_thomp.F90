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

module mod_cloud_thomp

  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams

  implicit none

  private

  real(rkx) , parameter :: entr = 0.35_rkx

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
  subroutine thomp_cldfrac(p,t,rho,qv,qc,qs,qi,iland,gridkm,cldfra)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: p , t , rho
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: qv , qc , qi , qs
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: iland
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: cldfra
    real(rkx) , intent(in) :: gridkm
    real(rkx) :: rh_00l , rh_00o , rh_00 , rhi_max
    real(rkx) , dimension(jci1:jci2,ici1:ici2,1:kz):: qvsat
    integer(ik4) :: i , j , k , kk
    real(rkx) :: tk , tc , qvsi , qvsw , rhum
    real(rkx) , dimension(kz) :: qvs1d , cfr1d , t1d , p1d
    real(rkx) , dimension(kz) :: r1d , qc1d , qi1d , qs1d

    ! First cut scale-aware. Higher resolution should require closer to
    ! saturated grid box for higher cloud fraction.  Simple functions
    ! chosen based on mocko and cotton (1995) starting point and desire
    ! to get near 100% rh as grid spacing moves toward 1.0km, but higher
    ! rh over ocean required as compared to over land.

    rh_00l = 0.839_rkx + sqrt(d_one/(50.0_rkx+gridkm*gridkm*gridkm*0.5_rkx))
    rh_00o = 0.879_rkx + sqrt(d_one/(100.0_rkx+gridkm*gridkm))

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          cldfra(j,i,k) = d_zero
          rhi_max = d_zero
          if ( qc(j,i,k) + qi(j,i,k) > 1.e-4_rkx ) then
            cldfra(j,i,k) = d_one
            qvsat(j,i,k) = qv(j,i,k)
          else
            tk = t(j,i,k)
            tc = tk - tzero
            qvsw = rslf(p(j,i,k), tk)
            qvsi = rsif(p(j,i,k), tk)
            if ( tc >= -12.0_rkx ) then
              qvsat(j,i,k) = qvsw
            else if ( tc < -30.0_rkx )  then
              qvsat(j,i,k) = qvsi
            else
              qvsat(j,i,k) = qvsw - &
                (qvsw-qvsi)*(-12.0_rkx-tc)/(-12.0_rkx+30.0_rkx)
            end if
            rhum = max(0.0_rkx,min(1.0_rkx,qv(j,i,k)/qvsat(j,i,k)))
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
              cldfra(j,i,k) = max(d_zero, ((rh_00-rhum)/(rh_00-rhi_max)) * &
                                          ((rh_00-rhum)/(rh_00-rhi_max)))
            end if
            cldfra(j,i,k) = max(d_zero, min(cldfra(j,i,k), d_one))
          end if
        end do
      end do
    end do

    ! Prepare for a 1-d column to find various cloud layers.

    do i = ici1 , ici2
      do j = jci1 , jci2
        do k = 1 , kz
          kk = kzp1-k
          qvs1d(kk) = qvsat(j,i,k)
          cfr1d(kk) = cldfra(j,i,k)
          t1d(kk) = t(j,i,k)
          p1d(kk) = p(j,i,k)
          r1d(kk) = rho(j,i,k)
          qc1d(kk) = qc(j,i,k)
          qi1d(kk) = qi(j,i,k)
          qs1d(kk) = qs(j,i,k)
        end do
        call find_cloudlayers(qvs1d,cfr1d,t1d,p1d,r1d,qc1d,qi1d,qs1d)
        do k = 1 , kz
          kk = kzp1-k
          cldfra(j,i,k) = cfr1d(kk)
          !qc(j,i,k) = qc1d(kk)
          !qi(j,i,k) = qi1d(kk)
        end do
      end do
    end do

    contains
     !
     ! This function calculates the liquid saturation vapor mixing ratio as
     ! a function of temperature and pressure
     !
     pure real(rkx) function rslf(p,t)
       implicit none
       real(rkx) , intent(in) :: p , t
       real(rkx) :: esl , x
       real(rkx) , parameter :: c0 =  0.611583699e03_rkx
       real(rkx) , parameter :: c1 =  0.444606896e02_rkx
       real(rkx) , parameter :: c2 =  0.143177157e01_rkx
       real(rkx) , parameter :: c3 =  0.264224321e-1_rkx
       real(rkx) , parameter :: c4 =  0.299291081e-3_rkx
       real(rkx) , parameter :: c5 =  0.203154182e-5_rkx
       real(rkx) , parameter :: c6 =  0.702620698e-8_rkx
       real(rkx) , parameter :: c7 =  0.379534310e-11_rkx
       real(rkx) , parameter :: c8 = -0.321582393e-13_rkx
       x = max(-80.0_rkx,t-tzero)
       esl = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
       rslf = ep2 * esl/(p-esl)
     end function rslf
     !
     ! This function calculates the ice saturation vapor mixing ratio as a
     ! function of temperature and pressure
     !
     pure real(rkx) function rsif(p,t)
       implicit none
       real(rkx) , intent(in) :: p , t
       real(rkx) :: esi , x
       real(rkx) , parameter :: c0 = 0.609868993e03_rkx
       real(rkx) , parameter :: c1 = 0.499320233e02_rkx
       real(rkx) , parameter :: c2 = 0.184672631e01_rkx
       real(rkx) , parameter :: c3 = 0.402737184e-1_rkx
       real(rkx) , parameter :: c4 = 0.565392987e-3_rkx
       real(rkx) , parameter :: c5 = 0.521693933e-5_rkx
       real(rkx) , parameter :: c6 = 0.307839583e-7_rkx
       real(rkx) , parameter :: c7 = 0.105785160e-9_rkx
       real(rkx) , parameter :: c8 = 0.161444444e-12_rkx
       x = max(-80.0_rkx,t-tzero)
       esi = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
       rsif = ep2 * esi/(p-esi)
     end function rsif

  end subroutine thomp_cldfrac
  !
  ! From cloud fraction array, find clouds of multi-level depth and compute
  ! a reasonable value of lwp or iwp that might be contained in that depth,
  ! unless existing lwc/iwc is already there.
  !
  subroutine find_cloudlayers(qvs1d,cfr1d,t1d,p1d,r1d,qc1d,qi1d,qs1d)
    implicit none
    real(rkx) , dimension(kz) , intent(in) :: qvs1d , t1d , p1d , r1d
    real(rkx) , dimension(kz) , intent(inout) :: cfr1d
    real(rkx) , dimension(kz) , intent(inout) :: qc1d , qi1d , qs1d

    real(rkx) , dimension(kz) :: theta , dz
    real(rkx) :: z1 , z2 , theta1 , theta2 , ht1 , ht2
    integer(ik4) :: k , k2 , k_tropo , k_m12c , k_m40c , k_cldb , k_cldt , kbot
    logical :: in_cloud

    k_m12c = 0
    k_m40c = 0
    do k = kz , 1 , -1
      theta(k) = t1d(k)*((100000.0_rkx/p1d(k))**rovcp)
      if ( t1d(k)-tzero > -40.0_rkx ) k_m40c = max(k_m40c, k)
      if ( t1d(k)-tzero > -12.0_rkx ) k_m12c = max(k_m12c, k)
    end do
    if ( k_m40c <= 1 ) k_m40c = 1
    if ( k_m12c <= 1 ) k_m12c = 1

    z2 = 44307.692_rkx * (d_one - (p1d(kz)/stdp)**0.190_rkx)
    do k = kz-1 , 1 , -1
      z1 = 44307.692_rkx * (d_one - (p1d(k)/stdp)**0.190_rkx)
      dz(k+1) = z2 - z1
      z2 = z1
    end do
    dz(1) = dz(2)
    !
    ! Find tropopause height, best surrogate, because we would not really
    ! wish to put fake clouds into the stratosphere.  The 10/1500 ratio
    ! d(theta)/d(z) approximates a vertical line on typical skewt chart
    ! near typical (mid-latitude) tropopause height.  Since messy data
    ! could give us a false signal of such a transition, do the check over
    ! three k-level change, not just a level-to-level check.  This method
    ! has potential failure in arctic-like conditions with extremely low
    ! tropopause height, as would any other diagnostic, so ensure resulting
    ! k_tropo level is above 4km.
    !
    do k = kz-3 , 1 , -1
      theta1 = theta(k)
      theta2 = theta(k+2)
      ht1 = 44307.692_rkx * (d_one - (p1d(k)/stdp)**0.190_rkx)
      ht2 = 44307.692_rkx * (d_one - (p1d(k+2)/stdp)**0.190_rkx)
      if ( (((theta2-theta1)/(ht2-ht1)) < 10.0_rkx/1500.0_rkx ) .and. &
           (ht1 < 19000.0_rkx) .and. (ht1 > 4000.0_rkx) ) then
        exit
      end if
    end do
    k_tropo = max(3, k+2)
    !
    ! Eliminate possible fractional clouds above supposed tropopause.
    !
    do k = k_tropo+1 , kz
      if (cfr1d(k) > 0.0_rkx .and. cfr1d(k) < 0.999_rkx ) then
        cfr1d(k) = d_zero
      end if
    end do
    !
    ! We would like to prevent fractional clouds below lcl in idealized
    ! situation with deep well-mixed convective pbl, that otherwise is
    ! likely to get clouds in more realistic capping inversion layer.
    !
    kbot = 3
    do k = kbot , k_m12c
      if ( (theta(k)-theta(k-1)) > 0.05e-3_rkx*dz(k) ) exit
    end do
    kbot = max(2, k-1)
    do k = 1 , kbot
      if ( cfr1d(k) > 0.0_rkx .and. cfr1d(k) < 0.999_rkx ) cfr1d(k) = d_zero
    end do
    !
    ! Starting below tropo height, if cloud fraction greater than 1 percent,
    ! compute an approximate total layer depth of cloud, determine a total
    ! liquid water/ice path (lwp/iwp), then reduce that amount with tuning
    ! parameter to represent entrainment factor, then divide up lwp/iwp
    ! into delta-z weighted amounts for individual levels per cloud layer.
    !
    k_cldb = k_tropo
    in_cloud = .false.
    k = k_tropo
    do while ( .not. in_cloud .and. k > k_m12c )
      k_cldt = 0
      if ( cfr1d(k) >= lowcld ) then
        in_cloud = .true.
        k_cldt = max(k_cldt, k)
      end if
      if ( in_cloud ) then
        do k2 = k_cldt-1 , k_m12c , -1
          if ( cfr1d(k2) < lowcld .or. k2 == k_m12c ) then
            k_cldb = k2+1
            exit
          end if
        end do
        in_cloud = .false.
      end if
      if ( (k_cldt - k_cldb + 1) >= 2 ) then
        call adjust_cloudice(cfr1d,qi1d,qs1d,qvs1d,t1d,r1d,dz,k_cldb,k_cldt)
        k = k_cldb
      end if
      k = k - 1
    end do

    k_cldb = k_tropo
    in_cloud = .false.
    k = k_m12c
    do while ( .not. in_cloud .and. k > kbot )
      k_cldt = 0
      if ( cfr1d(k) >= lowcld ) then
        in_cloud = .true.
        k_cldt = max(k_cldt, k)
      end if
      if ( in_cloud ) then
        do k2 = k_cldt-1 , kbot , -1
          if ( cfr1d(k2) < lowcld .or. k2 == kbot) then
            k_cldb = k2+1
            exit
          end if
        end do
        in_cloud = .false.
      end if
      if ( (k_cldt - k_cldb + 1) >= 2 ) then
        call adjust_cloudh2o(cfr1d,qc1d,qvs1d,t1d,r1d,dz,k_cldb,k_cldt)
        k = k_cldb
      end if
      k = k - 1
    end do

    ! Do a final total column adjustment since we may have added more than 1mm
    ! lwp/iwp for multiple cloud decks.

    call adjust_cloudfinal(cfr1d,qc1d,qi1d,r1d,dz,k_tropo)

  end subroutine find_cloudlayers

  subroutine adjust_cloudice(cfr,qi,qs,qvs,t,rho,dz,k1,k2)
    implicit none
    integer(ik4) , intent(in):: k1 , k2
    real(rkx) , dimension(kz) , intent(in):: cfr , qvs , t , rho , dz
    real(rkx) , dimension(kz) , intent(inout):: qi , qs
    real(rkx) :: iwp , max_iwp , tdz , this_iwp , iwp_exists
    integer(ik4) :: k

    max_iwp = abs(qvs(k2-1)-qvs(k2))*rho(k2-1)*dz(k2-1)
    tdz = d_zero
    iwp = d_zero
    iwp_exists = d_zero
    do k = k1 , k2
      tdz = tdz + dz(k)
      iwp = iwp + max(d_zero, (qvs(k-1)-qvs(k))*rho(k)*dz(k))
      iwp_exists = iwp_exists + (qi(k)+qs(k))*rho(k)*dz(k)
    end do
    if ( iwp_exists > d_one ) return
    max_iwp = max(max_iwp*(d_one-entr), min(d_one, iwp*(d_one-entr)))
    do k = k1 , k2
      this_iwp = max_iwp*dz(k)/tdz
      if (cfr(k) > lowcld .and. cfr(k) < hicld .and. t(k) >= 203.16_rkx ) then
        qi(k) = qi(k) + cfr(k)*cfr(k)*this_iwp/rho(k)/dz(k)
      end if
    end do
  end subroutine adjust_cloudice

  subroutine adjust_cloudh2o(cfr,qc,qvs,t,rho,dz,k1,k2)
    implicit none
    integer(ik4) , intent(in):: k1 , k2
    real(rkx) , dimension(kz) :: cfr , qc , qvs , t , rho , dz
    real(rkx) :: lwp , max_lwp , tdz , this_lwp , lwp_exists
    integer(ik4) :: k

    max_lwp = abs(qvs(k2-1)-qvs(k2))*rho(k2-1)*dz(k2-1)
    tdz = d_zero
    lwp = d_zero
    lwp_exists = d_zero
    do k = k1 , k2
      tdz = tdz + dz(k)
      lwp = lwp + max(d_zero, (qvs(k-1)-qvs(k))*rho(k)*dz(k))
      lwp_exists = lwp_exists + qc(k)*rho(k)*dz(k)
    end do
    if ( lwp_exists > d_one ) return
    max_lwp = max(max_lwp*(d_one-entr), min(d_one, lwp*(d_one-entr)))
    do k = k1 , k2
      this_lwp = max_lwp*dz(k)/tdz
      if ( cfr(k) > 0.95_rkx .and. &
           qc(k) < 1.e-7_rkx .and. t(k) < 253.16_rkx ) then
        qc(k) = qc(k) + 0.05_rkx*this_lwp/rho(k)/dz(k)
      else if ( cfr(k) > lowcld .and. cfr(k) < hicld .and. &
                t(k) < tzero .and. t(k) >= 253.16_rkx ) then
        qc(k) = qc(k) + cfr(k)*this_lwp/rho(k)/dz(k)
      else if ( cfr(k) > lowcld .and. cfr(k) < hicld .and. &
                t(k) <= 298.16_rkx .and.t(k) >= tzero ) then
        qc(k) = qc(k) + cfr(k)*cfr(k)*this_lwp/rho(k)/dz(k)
      end if
    end do

  end subroutine adjust_cloudh2o
  !
  ! Do not alter any grid-explicitly resolved hydrometeors, rather only
  ! the supposed amounts due to the cloud fraction scheme.
  !
  subroutine adjust_cloudfinal(cfr,qc,qi,rho,dz,k_tropo)
    implicit none
    integer(ik4) , intent(in) :: k_tropo
    real(rkx) , dimension(kz) , intent(in) :: cfr , rho , dz
    real(rkx) , dimension(kz) , intent(inout) :: qc , qi
    real(rkx) :: lwp , iwp , xfac
    integer(ik4) :: k

    lwp = d_zero
    do k = 1 , k_tropo
      if ( cfr(k) > lowcld .and. cfr(k) < hicld ) then
        lwp = lwp + qc(k)*rho(k)*dz(k)
      end if
    end do
    iwp = d_zero
    do k = 1 , k_tropo
      if ( cfr(k) > lowcld .and. cfr(k) < hicld ) then
        iwp = iwp + qi(k)*rho(k)*dz(k)
      end if
    end do
    if ( lwp > d_one ) then
      xfac = d_one/lwp
      do k = 1 , k_tropo
        if ( cfr(k) > lowcld .and. cfr(k) < hicld ) then
          qc(k) = qc(k)*xfac
        end if
      end do
    end if
    if ( iwp > d_one ) then
      xfac = d_one/iwp
      do k = 1 , k_tropo
        if ( cfr(k) > lowcld .and. cfr(k) < hicld ) then
          qi(k) = qi(k)*xfac
        end if
      end do
    end if
  end subroutine adjust_cloudfinal

end module mod_cloud_thomp
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
