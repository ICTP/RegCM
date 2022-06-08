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

module mod_cu_kuo
  !
  ! This module implements Kuo cumulus parameterization scheme.
  ! The basic method follows Anthes and Keyser (1979) and Kuo (1983).
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_mppparam
  use mod_cu_common
  use mod_constants
  use mod_service
  use mod_runparams , only : iqv , dt , ichem , dsigma , hsigma , qcon
  use mod_regcm_types

  implicit none

  private

  public :: allocate_mod_cu_kuo , cupara
  !
  ! qdcrit : the precipitation threshold for moisture convergence.
  ! pert   : perturbation temperature
  ! perq   : perturbation mixing ratio
  ! dlt    : temperature difference used to allow over shooting.
  ! cdscld : critical cloud depth in delta sigma.
  !
  real(rkx) , parameter :: qdcrit = 3.0e-7_rkx
  real(rkx) , parameter :: pert   = 1.0_rkx
  real(rkx) , parameter :: perq   = 1.0e-3_rkx
  real(rkx) , parameter :: dlt    = 3.0_rkx
  real(rkx) , parameter :: cdscld = 0.3_rkx
  real(rkx) , parameter :: bfac   = 0.5_rkx

  real(rkx) , public , pointer , dimension(:) :: qwght
  real(rkx) , public , pointer , dimension(:,:,:) :: twght , vqflx

  integer(ik4) , public :: k700

  real(rkx) , parameter :: svpt0 = tzero
  real(rkx) , parameter :: svp1 = 0.6112_rkx
  real(rkx) , parameter :: svp3 = 29.65_rkx
  real(rkx) , parameter :: svp2 = 17.67_rkx
  real(rkx) , parameter :: tauht = 7200.0_rkx

  contains

  subroutine allocate_mod_cu_kuo
    implicit none
    call getmem1d(qwght,1,kz,'cu_kuo:qwght')
    call getmem3d(twght,1,kz,5,kz,1,kz-3,'cu_kuo:twght')
    call getmem3d(vqflx,1,kz,5,kz,1,kz-3,'cu_kuo:vqflx')
  end subroutine allocate_mod_cu_kuo

  subroutine cupara(m2c)
    implicit none
    type(mod_2_cum) , intent(in) :: m2c
    real(rkx) :: apcnt , arh , c301 , dalr , deqt , dlnp , dplr , dsc ,   &
            ee , eddyf , emax , eqt , eqtm , plcl , pmax , pratec , psg , &
            psx , q , qmax , qs , rh , sca , siglcl , ff ,  suma , sumb , &
            t1 , tdmax , tdpt , tlcl , tmax , tmean , ttconv , ttp ,      &
            ttsum , xsav , zlcl
    integer(ik4) :: i , j , k , kbase , kk , ktop
    logical :: lconv
    real(rkx) , dimension(kz) :: tux , pux , qux , seqt
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cupara'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    !
    ! compute the moisture convergence in a column:
    ! at this stage, qten(j,i,k,iqv) only includes horizontal advection.
    ! sca: is the amount of total moisture convergence
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        lconv = .false.
        psx = m2c%psf(j,i) * d_r1000 ! Put in cb
        pux(:) = m2c%pas(j,i,:) * d_r1000 ! Put in cb
        tux(:) = m2c%tas(j,i,:)
        qux(:) = m2c%qxas(j,i,:,iqv)
        sca = d_zero
        do k = 1 , kz
          sca = sca + m2c%dynqx(j,i,k,iqv) * dsigma(k)
        end do
        !
        ! determine if moist convection exists:
        !
        if ( sca >= qdcrit ) then
          !
          ! check for stability
          !
          ! 1) compute eqt (equivalent potential temperature)
          !    between surface and 700 mb, with perturbation temperature
          !    and moisture added. the maximum eqt will be regarded
          !    as the origin of air parcel that produce cloud.
          !
          eqtm = d_zero
          pmax = pux(k700)
          qmax = qux(k700)
          tmax = tux(k700)
          do k = k700 , kz
            psg = pux(k)
            ttp = tux(k) + pert
            q = qux(k) + perq
            t1 = ttp*(d_100/psg)**rovcp
            ee = psg * q/(0.622_rkx+q)
            tdpt = min((d_one/(d_one/svpt0-rwat/wlhv*log(ee/0.611_rkx))),ttp)
            tlcl = tdpt - (0.212_rkx + 1.571e-3_rkx*(tdpt-svpt0) - &
                           4.36e-4_rkx*(ttp-svpt0)) * (ttp-tdpt)
            eqt = t1 * exp(wlhvocp*q/tlcl)
            if ( eqt > eqtm ) then
              eqtm = eqt
              tmax = ttp
              qmax = q
              pmax = psg
            end if
          end do
          !
          ! 2) compute lcl, get the sigma and p of lcl
          !
          emax = qmax*pmax/(ep2+qmax)
          tdmax = (svp3*log(emax/svp1)-svp2*svpt0)/(log(emax/svp1)-svp2)
          dalr = egrav*rcpd
          dplr = (egrav*tdmax*tdmax)/(ep2*wlhv*tmax)
          zlcl = (tmax-tdmax)/(dalr-dplr)
          tlcl = tmax - dalr*zlcl
          tmean = d_half*(tmax+tlcl)
          dlnp = (egrav*zlcl)/(rgas*tmean)
          plcl = pmax*exp(-dlnp)
          siglcl = (plcl-ptop)/psx
          !
          ! 3) compute seqt (saturation equivalent potential temperature)
          !    of all the levels that are above the lcl
          !
          do k = kz - 1 , 1 , -1
            kbase = k
            if ( hsigma(k) <= siglcl ) exit
          end do
          !
          ! kbase is the layer where lcl is located.
          !
          seqt(:) = d_zero
          do k = 1 , kbase
            ttp = tux(k)
            psg = pux(k)
            t1 = ttp*(d_100/psg)**rovcp
            ee = svp1*exp(svp2*(ttp-svpt0)/(ttp-svp3))
            qs = ep2*ee/(psg-ee)
            seqt(k) = t1*exp(wlhvocp*qs/ttp)
          end do
          !
          ! 4) when seqt = eqt + deqt, cloud top is reached.
          !    eqt is the eqt of cloud (same as lcl eqt).
          !
          ktop = max(kbase-3,1)
          do k = 1 , kbase
            kk = kbase + 1 - k
            deqt = seqt(kk) - eqtm
            if ( deqt > dlt ) exit
          end do
          !
          ! cloud top has been reached
          !
          ktop = min(kk,ktop)
          !
          ! 5) check cloud depth
          !    if cloud depth is less than critical depth (cdscld = 0.3),
          !    the convection is killed
          !
          dsc = (siglcl-hsigma(ktop))
          if ( dsc >= cdscld ) then
            !
            ! 6) check negative area
            !    if negative area is larger than the positive area
            !    convection is killed.
            !
            ttsum = d_zero
            do k = ktop , kbase
              ttsum = ttsum + (eqtm-seqt(k))*dsigma(k)
            end do

            if ( ttsum >= d_zero ) then
              !
              ! convection exist, compute convective flux of water vapor and
              ! latent heating
              ! icon   : is a counter which keep track the total points
              !          where deep convection occurs.
              ! c301   : is the 'b' factor in kuo's scheme.
              !
              lconv = .true.
              total_precip_points = total_precip_points + 1
              suma = d_zero
              sumb = d_zero
              arh = d_zero
              qwght(:) = d_zero
              do k = ktop , kz
                psg = pux(k)
                ttp = tux(k)
                ee = svp1*exp(svp2*(tux(k)-svpt0)/(tux(k)-svp3))
                qs = ep2*ee/(pux(k)-ee)
                rh = qux(k)/qs
                rh = max(min(rh,d_one),d_zero)
                xsav = (d_one-rh)*qs
                qwght(k) = xsav
                sumb = sumb + qs*dsigma(k)
                suma = suma + xsav*dsigma(k)
                arh = arh + rh*qs*dsigma(k)
              end do
              arh = arh/sumb
              c301 = bfac*(d_one-arh)
              if ( c301 < d_zero ) c301 = d_zero
              if ( c301 > d_one ) c301 = d_one
              if ( suma <= d_zero ) then
                c301 = d_zero
                suma = d_one
              end if
              do k = ktop , kz
                qwght(k) = qwght(k)/suma
              end do
              do k = ktop , kbase
                ttconv = wlhvocp*(d_one-c301)*twght(k,kbase,ktop)*sca
                apcnt = (d_one-c301)*sca/4.3e-3_rkx
                eddyf = apcnt*vqflx(k,kbase,ktop)
                cu_qten(j,i,k,iqv) = (c301*qwght(k)*sca + eddyf) / psx
                cu_tten(j,i,k) = ttconv / psx
              end do
              cu_ktop(j,i) = ktop
              cu_kbot(j,i) = kbase
              ! the unit for rainfall is kg m-2 s-1
              pratec = (d_one-c301)*sca*d_100*regrav
              if ( pratec > dlowval ) then
                ! instantaneous precipitation rate for use in surface (mm/s)
                cu_prate(j,i) = cu_prate(j,i) + pratec
                if ( ichem == 1 ) then
                  ! build for chemistry 3d table of cons precipitation rate
                  ! from the surface to the top of the convection
                  do k = 1 , ktop-1
                    cu_convpr(j,i,k) = pratec
                  end do
                end if
              end if
            end if
          end if
        end if
        if ( .not. lconv ) then
          !
          ! convection do not exist, compute the vertical advection term:
          !
          do k = 2 , kz
            if ( m2c%qq1(j,i,k)   / m2c%psa(j,i) > minqq .and. &
                 m2c%qq1(j,i,k-1) / m2c%psa(j,i) > minqq ) then
              ff = ((m2c%qq1(j,i,k)*(m2c%qq1(j,i,k-1) / &
                m2c%qq1(j,i,k))**qcon(k)) * m2c%qdot(j,i,k))/m2c%psa(j,i)
              cu_qten(j,i,k-1,iqv) = cu_qten(j,i,k-1,iqv) - ff/dsigma(k-1)
              cu_qten(j,i,k,iqv)   = cu_qten(j,i,k,iqv) + ff/dsigma(k)
            end if
          end do
        end if
      end do
    end do

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cupara

end module mod_cu_kuo

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
