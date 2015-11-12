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
  use mod_service
  use mod_runparams , only : iqv , dt , dtsec , ichem , dsigma , hsigma , qcon
  use mod_regcm_types

  implicit none

  private

  public :: allocate_mod_cu_kuo , cupara , htdiff
  !
  ! qdcrit : the precipitation threshold for moisture convergence.
  ! pert   : perturbation temperature
  ! perq   : perturbation mixing ratio
  ! dlt    : temperature difference used to allow over shooting.
  ! cdscld : critical cloud depth in delta sigma.
  !
  real(rk8) , parameter :: qdcrit = 3.0D-7
  real(rk8) , parameter :: pert   = 1.0D0
  real(rk8) , parameter :: perq   = 1.0D-3
  real(rk8) , parameter :: dlt    = 3.0D0
  real(rk8) , parameter :: cdscld = 0.3D0

  real(rk8) , public , pointer , dimension(:,:,:) :: rsheat , rswat
  real(rk8) , public , pointer , dimension(:) :: qwght
  real(rk8) , public , pointer , dimension(:,:,:) :: twght , vqflx
  real(rk8) , pointer , dimension(:,:,:) :: wrkkuo1
  real(rk8) , pointer , dimension(:,:,:) :: wrkkuo2

  integer(ik4) , public :: k700

  contains

  subroutine allocate_mod_cu_kuo
    implicit none
    call getmem3d(rsheat,jci1,jci2,ici1,ici2,1,kz,'cu_kuo:rsheat')
    call getmem3d(rswat,jci1,jci2,ici1,ici2,1,kz,'cu_kuo:rswat')
    call getmem1d(qwght,1,kz,'cu_kuo:qwght')
    call getmem3d(twght,1,kz,5,kz,1,kz-3,'cu_kuo:twght')
    call getmem3d(vqflx,1,kz,5,kz,1,kz-3,'cu_kuo:vqflx')
    call getmem3d(wrkkuo1,jce1-ma%jbl1,jce2+ma%jbr1, &
                          ice1-ma%ibb1,ice2+ma%ibt1,1,kz,'tendency:wrkkuo1')
    call getmem3d(wrkkuo2,jce1-ma%jbl1,jce2+ma%jbr1, &
                          ice1-ma%ibb1,ice2+ma%ibt1,1,kz,'tendency:wrkkuo2')
  end subroutine allocate_mod_cu_kuo

  subroutine cupara(m2c,c2m)
    implicit none
    type(mod_2_cum) , intent(in) :: m2c
    type(cum_2_mod) , intent(inout) :: c2m
    real(rk8) :: apcnt , arh , c301 , dalr , deqt , dlnp , dplr , dsc ,   &
            eddyf , emax , eqt , eqtm , plcl , pmax , pratec ,            &
            q , qmax , qs , rh , rsht , rswt , sca , siglcl ,             &
            suma , sumb , t1 , tdmax , tlcl , tmax , tmean , ttconv ,     &
            ttp , ttsum , xsav , zlcl
    integer(ik4) :: i , j , k , kbase , kbaseb , kk , ktop
    real(rk8) , dimension(kz) :: seqt
    real(rk8) , dimension(kz) :: tmp3
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cupara'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    !
    ! c2m%kcumtop = top level of cumulus clouds
    ! c2m%kcumbot = bottom level of cumulus clouds
    !
    c2m%kcumtop(:,:) = 0
    c2m%kcumbot(:,:) = 0
    if ( ichem == 1 ) c2m%convpr(:,:,:) = d_zero
    !
    ! compute the moisture convergence in a column:
    ! at this stage, c2m%qxten(j,i,k,iqv) only includes horizontal advection.
    ! sca: is the amount of total moisture convergence
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        sca = d_zero
        do k = 1 , kz
          sca = sca + c2m%qxten(j,i,k,iqv)*dsigma(k)
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
          pmax =  m2c%pas(j,i,kz)
          qmax = m2c%qxas(j,i,kz,iqv)
          tmax =  m2c%tas(j,i,kz)
          do k = k700 , kz
            ttp = m2c%tas(j,i,k) + pert
            q = m2c%qxas(j,i,k,iqv) + perq
            t1 = ttp*(1.0D5/m2c%pas(j,i,k))**rovcp
            eqt = t1*exp(wlhvocp*q/ttp)
            if ( eqt > eqtm ) then
              eqtm = eqt
              tmax = ttp
              qmax = q
              pmax = m2c%pas(j,i,k)*d_r100
            end if
          end do
          !
          ! 2) compute lcl, get the sigma and p of lcl
          !
          emax = qmax*pmax/(ep2+qmax)
          tdmax = 5418.12D0/(19.84659D0-log(emax/0.611D0))
          dalr = egrav*rcpd
          dplr = (egrav*tdmax*tdmax)/(ep2*wlhv*tmax)
          zlcl = (tmax-tdmax)/(dalr-dplr)
          tlcl = tmax - dalr*zlcl
          tmean = (tmax+tlcl)*d_half
          dlnp = (egrav*zlcl)/(rgas*tmean)
          plcl = pmax*exp(-dlnp)
          siglcl = (plcl-ptop)/m2c%psb(j,i)
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
          do k = 1 , kbase
            ttp = m2c%tas(j,i,k)
            qs = pfqsat(ttp,m2c%pas(j,i,k))
            t1 = ttp*(1.0D5/m2c%pas(j,i,k))**rovcp
            seqt(k) = t1*exp(wlhvocp*rcpd*qs/ttp)
          end do
          !
          ! 4) when seqt = eqt + dt, cloud top is reached.
          !    eqt is the eqt of cloud (same as lcl eqt).
          !
          do kk = 1 , kbase
            k = kbase + 1 - kk
            deqt = seqt(k) - eqtm
            if ( deqt > dlt ) exit
          end do
          !
          ! cloud top has been reached
          !
          ktop = min(kbase-3,k)
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
              ttsum = (eqtm-seqt(k))*dsigma(k) + ttsum
            end do
            if ( ttsum >= d_zero ) then
              !
              ! convection exist, compute convective flux of water vapor and
              ! latent heating
              ! icon   : is a counter which keep track the total points
              !          where deep convection occurs.
              ! c301   : is the 'b' factor in kuo's scheme.
              !
              total_precip_points = total_precip_points + 1
              suma = d_zero
              sumb = d_zero
              arh = d_zero
              do k = 1 , kz
                qwght(k) = d_zero
              end do
              do k = ktop , kz
                ttp = m2c%tas(j,i,k)
                qs = pfqsat(ttp,m2c%pas(j,i,k))
                rh = m2c%qxas(j,i,k,iqv)/qs
                rh = max(min(rh,d_one),d_zero)
                xsav = (d_one-rh)*qs
                qwght(k) = xsav
                sumb = sumb + qs*dsigma(k)
                arh = arh + rh*qs*dsigma(k)
                suma = suma + xsav*dsigma(k)
              end do
              arh = arh/sumb
              c301 = d_two*(d_one-arh)
              if ( c301 < d_zero ) c301 = d_zero
              if ( c301 > d_one ) c301 = d_one
              if ( suma <= d_zero ) then
                c301 = d_zero
                suma = d_one
              end if
              do k = ktop , kz
                qwght(k) = qwght(k)/suma
              end do
              do k = 1 , kz
                ttconv = wlhvocp*(d_one-c301)*twght(k,kbase,ktop)*sca
                rsheat(j,i,k) = rsheat(j,i,k) + ttconv*dt
                apcnt = (d_one-c301)*sca/4.3D-3
                eddyf = apcnt*vqflx(k,kbase,ktop)
                c2m%qxten(j,i,k,iqv) = eddyf
                rswat(j,i,k) = rswat(j,i,k) + c301*qwght(k)*sca*dt
              end do
              kbaseb = min0(kbase,kzm2)
              c2m%kcumtop(j,i) = ktop
              c2m%kcumbot(j,i) = kbaseb
              ! the unit for rainfall is mm.
              pratec = (d_one-c301)*sca*d_1000*regrav
              if ( pratec > dlowval ) then
                c2m%rainc(j,i) = c2m%rainc(j,i) + pratec*dtsec
                ! instantaneous precipitation rate for use in surface (mm/s)
                c2m%pcratec(j,i) = c2m%pcratec(j,i) + pratec
              end if
              if ( ichem == 1 ) then
                ! build for chemistry 3d table of cons precipitation rate
                ! from the surface to the top of the convection
                do k = 1 , ktop-1
                  c2m%convpr(j,i,kz-k+1) = pratec
                end do
              end if
              cycle
            end if
          end if
        end if
        !
        ! convection not exist, compute the vertical advection term:
        !
        tmp3(1) = d_zero
        do k = 2 , kz
          if ( m2c%qxas(j,i,k,iqv) < 1.0D-8 ) then
            tmp3(k) = d_zero
          else
            tmp3(k) = m2c%qxas(j,i,k,iqv) * &
                      (m2c%qxas(j,i,k-1,iqv)/m2c%qxas(j,i,k,iqv))**qcon(k)
          end if
        end do
        c2m%qxten(j,i,1,iqv) = c2m%qxten(j,i,1,iqv) - &
                m2c%qdot(j,i,2)*tmp3(2)/dsigma(1)
        do k = 2 , kzm1
          c2m%qxten(j,i,k,iqv) = c2m%qxten(j,i,k,iqv) - &
                  (m2c%qdot(j,i,k+1)*tmp3(k+1) - &
                   m2c%qdot(j,i,k)*tmp3(k))/dsigma(k)
        end do
        c2m%qxten(j,i,kz,iqv) = c2m%qxten(j,i,kz,iqv) + &
                m2c%qdot(j,i,kz)*tmp3(kz)/dsigma(kz)
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          rsheat(j,i,k) = max(rsheat(j,i,k),d_zero)
          rswat(j,i,k) = max(rswat(j,i,k),d_zero)
          rsht = rsheat(j,i,k)/tauht
          rswt = rswat(j,i,k)/tauht
          c2m%tten(j,i,k) = c2m%tten(j,i,k) + rsht
          c2m%qxten(j,i,k,iqv) = c2m%qxten(j,i,k,iqv) + rswt
          rsheat(j,i,k) = rsheat(j,i,k)*(d_one-dtsec/tauht)
          rswat(j,i,k) = rswat(j,i,k)*(d_one-dtsec/tauht)
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine cupara

  subroutine htdiff(dxsq,akht1)
    implicit none
    real(rk8) , intent(in) :: akht1 , dxsq
    integer(ik4) :: i , j , k

    wrkkuo1(jci1:jci2,ici1:ici2,:) = rsheat(:,:,:)
    wrkkuo2(jci1:jci2,ici1:ici2,:) = rswat(:,:,:)
    if ( ma%has_bdyleft ) then
      wrkkuo1(jce1,ici1:ici2,:) = wrkkuo1(jci1,ici1:ici2,:)
      wrkkuo2(jce1,ici1:ici2,:) = wrkkuo2(jci1,ici1:ici2,:)
    end if
    if ( ma%has_bdyright ) then
      wrkkuo1(jce2,ici1:ici2,:) = wrkkuo1(jci2,ici1:ici2,:)
      wrkkuo2(jce2,ici1:ici2,:) = wrkkuo2(jci2,ici1:ici2,:)
    end if
    if ( ma%has_bdybottom ) then
      wrkkuo1(jci1:jci2,ice1,:) = wrkkuo1(jci1:jci2,ici1,:)
      wrkkuo2(jci1:jci2,ice1,:) = wrkkuo2(jci1:jci2,ici1,:)
    end if
    if ( ma%has_bdytop ) then
      wrkkuo1(jci1:jci2,ice2,:) = wrkkuo1(jci1:jci2,ici2,:)
      wrkkuo2(jci1:jci2,ice2,:) = wrkkuo2(jci1:jci2,ici2,:)
    end if
    if ( ma%has_bdyleft .and. ma%has_bdybottom ) then
      wrkkuo1(jce1,ice1,:) = wrkkuo1(jci1,ici1,:)
      wrkkuo2(jce1,ice1,:) = wrkkuo2(jci1,ici1,:)
    end if
    if ( ma%has_bdyleft .and. ma%has_bdytop ) then
      wrkkuo1(jce1,ice2,:) = wrkkuo1(jci1,ici2,:)
      wrkkuo2(jce1,ice2,:) = wrkkuo2(jci1,ici2,:)
    end if
    if ( ma%has_bdyright .and. ma%has_bdybottom ) then
      wrkkuo1(jce2,ice1,:) = wrkkuo1(jci2,ici1,:)
      wrkkuo2(jce2,ice1,:) = wrkkuo2(jci2,ici1,:)
    end if
    if ( ma%has_bdyright .and. ma%has_bdytop ) then
      wrkkuo1(jce2,ice2,:) = wrkkuo1(jci2,ici2,:)
      wrkkuo2(jce2,ice2,:) = wrkkuo2(jci2,ici2,:)
    end if
    call exchange(wrkkuo1,1,jce1,jce2,ice1,ice2,1,kz)
    call exchange(wrkkuo2,1,jce1,jce2,ice1,ice2,1,kz)
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          rsheat(j,i,k) = rsheat(j,i,k)+akht1*dt/dxsq * &
                   (wrkkuo1(j,i-1,k)+wrkkuo1(j,i+1,k) + &
                    wrkkuo1(j-1,i,k)+wrkkuo1(j+1,i,k)-d_four*wrkkuo1(j,i,k))
        end do
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          rswat(j,i,k) = rswat(j,i,k)+akht1*dt/dxsq * &
                (wrkkuo2(j,i-1,k)+wrkkuo2(j,i+1,k) + &
                 wrkkuo2(j-1,i,k)+wrkkuo2(j+1,i,k)-d_four*wrkkuo2(j,i,k))
        end do
      end do
    end do
  end subroutine htdiff
!
end module mod_cu_kuo
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
