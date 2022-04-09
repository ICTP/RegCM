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

module mod_micro_subex
  !
  ! Large Scale Precipitation computation
  ! Fractional cloud coverage and liquid water content calculation
  ! Heating term for explicit moisture scheme
  !
  ! -- Pal et al. 2000 JGR-Atmos
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams
  use mod_memutil
  use mod_mpmessage
  use mod_regcm_types

  implicit none

  private
  !
  ! Precip sum beginning from top
  !
  real(rkx) , public , pointer , dimension(:,:) :: qck1 , cgul , &
    cevap , xcevap , caccr

  real(rkx) :: maxlat
  real(rkx) , pointer , dimension(:,:) :: pptsum
  real(rkx) , pointer , dimension(:,:,:) :: dqc

  logical :: l_lat_hack = .false.
  public :: allocate_subex , init_subex , subex

  real(rkx) , parameter :: remfrc = 0.0_rkx
  real(rkx) , parameter :: rhow = 1000.0_rkx
  real(rkx) , parameter :: pptmin = 0.0_rkx
  real(rkx) , parameter :: actcld = 0.01_rkx
  real(rkx) , parameter :: actliq = 1.0e-8_rkx
  real(rkx) , parameter :: accrfrc = 0.5_rkx

  contains

  subroutine allocate_subex
    implicit none
    ! Those not. Note the external, internal change.
    call getmem2d(qck1,jci1,jci2,ici1,ici2,'subex:qck1')
    call getmem2d(cgul,jci1,jci2,ici1,ici2,'subex:cgul')
    call getmem2d(cevap,jci1,jci2,ici1,ici2,'subex:cevap')
    if ( l_lat_hack ) then
      call getmem2d(xcevap,jci1,jci2,ici1,ici2,'subex:xcevap')
    else
      call assignpnt(cevap,xcevap)
    end if
    call getmem2d(caccr,jci1,jci2,ici1,ici2,'subex:caccr')
    call getmem2d(pptsum,jci1,jci2,ici1,ici2,'subex:pptsum')
    call getmem3d(dqc,jci1,jci2,ici1,ici2,1,kz,'subex:dqc')
  end subroutine allocate_subex

  subroutine init_subex(xlat)
    use mod_mppparam , only : maxall
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat
    call maxall(maxval(xlat),maxlat)
  end subroutine init_subex
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                 c
  ! This subroutine computes the 'large scale' precipitation        c
  ! based on the excess of cloud water above a threshold value.     c
  ! The threshold (qcth) is temperature dependant and is based      c
  ! on in cloud measurements of liquid cloud water (not ice).       c
  ! Rain is only produced from the cloudy fraction of a grid cell   c
  ! but the calculated precip value is a grid cell average.         c
  !                                                                 c
  ! This routine also computes raindrop evaporation and accretion.  c
  !                                                                 c
  ! See Pal et al. 2000 JGR-Atmos for more information.             c
  !                                                                 c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine subex(mo2mc,mc2mo)
    implicit none
    type(mod_2_micro) , intent(in) :: mo2mc
    type(micro_2_mod) , intent(out) :: mc2mo
    real(rkx) :: dpovg , afc , pptacc , pptkm1 , pptmax ,       &
                pptnew , qcleft , qcw , qs , rdevap , qcincl ,  &
                rhcs , prainx , qcth , dqv , rlv , ocpm
    real(rkx) :: tcel
    integer(ik4) :: i , j , k , kk
    logical :: lsecind
    !
    ! 0. Compute dqc
    !
    lsecind = (ichem == 1 .and. iaerosol == 1 .and. iindirect == 2)
    if ( l_lat_hack ) then
      call sun_cevap
    end if
    if ( lsecind ) then
      if ( idiag > 0 ) then
        mc2mo%dia_qcr(:,:,:) = d_zero
        mc2mo%dia_qcl(:,:,:) = d_zero
      end if
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            afc = mc2mo%fcc(j,i,k)
            qcw = mo2mc%qcn(j,i,k)
            if ( qcw > actliq .and. afc > actcld ) then
              ! include aerosol second indirect effect on threshold
              ! auto-conversion
              ! rcrit is a critical cloud radius for cloud
              ! water undergoing autoconversion
              ! ccn = number of ccn /m3
              ! In cloud mixing ratio [kg/kg]
              qcincl = mo2mc%qcn(j,i,k)/afc
              qcth = mo2mc%ccn(j,i,k)*(4.0_rkx/3.0_rkx)*mathpi * &
                    ((rcrit*1e-6_rkx)**3)*rhow
              if ( idiag > 0 ) then
                mc2mo%dia_qcr(j,i,k) = qcth
                mc2mo%dia_qcl(j,i,k) = qcincl
              end if
              dqc(j,i,k) = max(qcincl - qcth,d_zero)
            else
              dqc(j,i,k) = d_zero
            end if
          end do
        end do
      end do
    else
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! 1ad. Implement here the formula for qcth.
            !   - Gultepe & Isaac, J. Clim, 1997, v10 p446 table 4, eq 5
            !   - The factor of 1000 converts from g/kg to kg/kg
            !   - The factor of cgul accounts for the fact that the Gultepe
            !     and Isaac equation is for mean cloud water while qcth is the
            !     theshhold for auto-conversion.
            afc = mc2mo%fcc(j,i,k)
            qcw = mo2mc%qcn(j,i,k)
            if ( qcw > actliq .and. afc > actcld ) then
              ! In cloud mixing ratio [kg/kg]
              qcincl = mo2mc%qcn(j,i,k)/afc
              tcel = mo2mc%t(j,i,k) - tzero   ![C][avg]
              qcth = cgul(j,i) * &
                 (d_10**(-0.48911_rkx+0.01344_rkx*tcel))*d_r1000
              ! Use same function of Lemus et al., 1997 as in lwc computation
              !qcth = cgul(j,i) * &
              !      clwfromt(mo2mc%t(j,i,k))/mo2mc%rho(j,i,k)*d_r1000
              dqc(j,i,k) = max(qcincl - qcth,d_zero)
            else
              dqc(j,i,k) = d_zero
            end if
          end do
        end do
      end do
    end if
    !--------------------------------------------------------------------
    ! 1. Compute the precipitation formed in each layer.
    !    The computations are performed from the top to the surface.
    !    - Auto-conversion: similar to SIMEX (Giorgi and Shields 1999).
    !    - Raindrop Accretion:  Beheng (1994); EQN 6
    !    - Raindrop Evaporation:  Sundqvist (1988); EQN 3.4a
    !--------------------------------------------------------------------
    !
    ! zero all accumulated precip
    !
    pptsum(:,:) = d_zero
    if ( ichem == 1 ) then
      mc2mo%remrat(:,:,:) = d_zero
      if ( lsecind .and. idiag > 0 ) then
        mc2mo%dia_acr(:,:,:) = d_zero
      end if
    end if
    ! 1a. Perform computations for the top layer (layer 1)
    !   maximum precipation rate (total cloud water/dt)
    do i = ici1 , ici2
      do j = jci1 , jci2
        pptnew = d_zero
        afc = min(mc2mo%fcc(j,i,1),d_one-actcld)                ![frac][avg]
        qcw = mo2mc%qcn(j,i,1)     ![kg/kg][avg]
        if ( qcw > actliq .and. afc > actcld ) then ! if there is a cloud
          pptmax = (d_one-remfrc)*qcw/dt    ![kg/kg/s][avg]
          ! 1ac. Compute the maximum precipation rate
          !      (i.e. total cloud water/dt) [kg/kg/s]
          ! 1ae. Compute the gridcell average autoconversion [kg/k g/s]
          pptnew = min(pptmax,qck1(j,i)*dqc(j,i,1)*afc)     ![kg/kg/s][avg]
          ! 1ag. Compute the amount of cloud water removed by raindrop
          !      accretion [kg/kg/s].  In the layer where the precipitation
          !      is formed, only half of the precipitation is assumed to
          !      accrete. 1aga. Compute the amount of water remaining in the
          !      cloud [kg/kg]
          qcleft = max(qcw - pptnew*dt,d_zero) ![kg/kg][avg]
          ! 1agb. Add fraction of the new precipitation can accrete.
          pptkm1 = accrfrc*pptnew/afc*mo2mc%rho(j,i,1)*dt  ![kg/m3][cld]
          ! 1agc. Accretion [kg/kg/s]=[m3/kg/s]*[kg/kg]*[kg/m3]
          pptacc = caccr(j,i)*qcleft*pptkm1  ![kg/kg/s][avg]
          ! 1agd. Update the precipitation accounting for the
          !       accretion [kg/kg/s]
          pptnew = min(pptmax,pptacc+pptnew) ![kg/kg/s][avg]
          ! 1ah. Accumulate precipitation and convert to kg/m2/s
          if ( pptnew > pptmin ) then
            dpovg = (mo2mc%pfs(j,i,2)-mo2mc%pfs(j,i,1))*regrav    ![kg/m2]
            pptsum(j,i) = pptnew*dpovg         ![kg/m2/s][avg]
            if ( idynamic == 3 ) then
              ! 1ai. Compute the cloud water tendency [kg/kg/s]
              ! [kg/kg/s][avg]
              mc2mo%qxten(j,i,1,iqc) = mc2mo%qxten(j,i,1,iqc) - pptnew
            else
              ! 1ai. Compute the cloud water tendency [kg/kg/s*cb]
              ! [kg/kg/s*cb][avg]
              mc2mo%qxten(j,i,1,iqc) = mc2mo%qxten(j,i,1,iqc) - &
                                    pptnew*mo2mc%psb(j,i)
            end if
            ! 1af. Compute the cloud removal rate (for chemistry) [1/s]
            if ( ichem == 1 ) then
              mc2mo%remrat(j,i,1) = pptnew/qcw
              if ( lsecind .and. idiag > 0 ) then
                mc2mo%dia_acr(j,i,1) = pptnew
              end if
            end if
          end if
        end if
      end do
    end do

    ! LAYER TWO TO KZ
    ! 1b. Perform computations for the 2nd layer to the surface
    ! precipation accumulated from above
    do k = 2 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          ! 1bb. Convert accumlated precipitation to kg/kg/s.
          !      Used for raindrop evaporation and accretion.
          dpovg = (mo2mc%pfs(j,i,k+1)-mo2mc%pfs(j,i,k))*regrav    ![kg/m2]
          afc = min(mc2mo%fcc(j,i,k),d_one-actcld)                ![frac][avg]
          qcw = mo2mc%qcn(j,i,k)                      ![kg/kg][avg]
          pptnew = d_zero
          if ( pptsum(j,i) > d_zero ) then
            pptkm1 = pptsum(j,i)/dpovg                ![kg/kg/s][avg]
          else
            pptkm1 = d_zero
          end if
          ! 1bc. Compute the raindrop evaporation in the clear portion of
          !      the gridcell.
          !  - It is assumed that raindrops do not evaporate in clouds
          !    and the rainfall from above is evenly distributed in
          !    gridcell (i.e. the gridcell average precipitation is used).
          if ( pptkm1 > pptmin ) then
            qs = pfwsat(mo2mc%t(j,i,k),mo2mc%phs(j,i,k))  ![kg/kg][avg]
            dqv = (qs-mo2mc%qvn(j,i,k))/dt
            if ( dqv > d_zero ) then
              ! 2bca. Compute the clear sky relative humidity
              rhcs = (mo2mc%rh(j,i,k)-afc*rhmax)/(d_one-afc)    ![frac][clr]
              rhcs = max(min(rhcs,rhmax),rhmin)                ![frac][clr]
              ! 2bcb. Raindrop evaporation [kg/kg/s]
              rdevap = xcevap(j,i)*(rhmax-rhcs)*sqrt(pptsum(j,i))*(d_one-afc)
              rdevap = min(rdevap,dqv)          ![kg/kg/s][avg]
              rdevap = min(rdevap,pptkm1)       ![kg/kg/s][avg]
              ! 2bcc. Update the precipitation accounting for the raindrop
              !       evaporation [kg/m2/s]
              if ( rdevap > dlowval ) then
                if ( rdevap*dpovg > pptsum(j,i) ) then
                  pptsum(j,i) = d_zero
                  rdevap = pptsum(j,i)/dpovg
                else
                  pptsum(j,i) = pptsum(j,i) - rdevap*dpovg ![kg/m2/s][avg]
                end if
                rlv = wlh(mo2mc%t(j,i,k))
                ocpm = d_one/(cpd*(d_one-mo2mc%qxx(j,i,k,iqv)) + &
                              cpv*mo2mc%qxx(j,i,k,iqv))
                pptkm1 = pptkm1 - rdevap
                if ( idynamic == 3 ) then
                  ! 2bcf. Compute the water vapor tendency [kg/kg/s]
                  ![kg/kg/s][avg]
                  mc2mo%qxten(j,i,k,iqv) = mc2mo%qxten(j,i,k,iqv) + rdevap
                  ! 2bcf. Compute the temperature tendency [K/s]
                  ![K/s][avg]
                  mc2mo%tten(j,i,k) = mc2mo%tten(j,i,k) - rlv*ocpm*rdevap
                else
                  ! 2bcf. Compute the water vapor tendency [kg/kg/s*cb]
                  ![kg/kg/s*cb][avg]
                  mc2mo%qxten(j,i,k,iqv) = mc2mo%qxten(j,i,k,iqv) + &
                                             rdevap*mo2mc%psb(j,i)
                  ! 2bcf. Compute the temperature tendency [K/s*cb]
                  ![k/s*cb][avg]
                  mc2mo%tten(j,i,k) = mc2mo%tten(j,i,k) - &
                                   rlv*ocpm*rdevap*mo2mc%psb(j,i)
                end if
              end if
            end if
          end if
          ! 1bd. Compute the autoconversion and accretion [kg/kg/s]
          if ( qcw > actliq .and. afc > actcld ) then ! if there is a cloud
            pptmax = (d_one-remfrc)*qcw/dt              ![kg/kg/s][avg]
            ! 1bdb. Compute the maximum precipation rate
            !       (i.e. total cloud water/dt) [kg/kg/s]
            ! 1bdd. Compute the gridcell average autoconversion [kg/kg/s]
            pptnew = min(pptmax,qck1(j,i)*dqc(j,i,k)*afc) ![kg/kg/s][avg]
            ! 1bf. Compute the amount of cloud water removed by raindrop
            !      accretion [kg/kg/s].  In the layer where the precipitation
            !      is formed, only half of the precipitation can accrete.
            ! 1bfa. Compute the amount of water remaining in the cloud [kg/kg]
            qcleft = max(qcw-pptnew*dt,d_zero)             ![kg/kg][avg]
            ! 1bfb. Add fraction of the new precipitation to the accumulated
            !       precipitation [kg/m3]
            ! [kg/m3][cld]
            pptkm1 = (pptkm1 + accrfrc*pptnew/afc) * mo2mc%rho(j,i,k) * dt
            ! 1bfc. accretion [kg/kg/s]
            pptacc = caccr(j,i)*qcleft*pptkm1              ![kg/kg/s][avg]
            ! 1bfd. Update the precipitation accounting for the
            !       accretion [kg/kg/s]
            pptnew = min(pptmax,pptacc+pptnew)
            if ( pptnew > pptmin ) then
              ! 1bg. Accumulate precipitation and convert to kg/m2/s
              pptsum(j,i) = pptsum(j,i) + pptnew*dpovg       ![kg/m2/s][avg]
              if ( idynamic == 3 ) then
                ! 1bh. Compute the cloud water tendency [kg/kg/s]
                ![kg/kg/s][avg]
                mc2mo%qxten(j,i,k,iqc) = mc2mo%qxten(j,i,k,iqc) - pptnew
              else
                ! 1bh. Compute the cloud water tendency [kg/kg/s*cb]
                ![kg/kg/s*cb][avg]
                mc2mo%qxten(j,i,k,iqc) = mc2mo%qxten(j,i,k,iqc) - &
                             pptnew*mo2mc%psb(j,i)
              end if
              ! 1be. Compute the cloud removal rate (for chemistry) [1/s]
              if ( ichem == 1 ) then
                mc2mo%remrat(j,i,k) = pptnew/qcw
                if ( lsecind .and. idiag > 0 ) then
                  mc2mo%dia_acr(j,i,k) = pptnew
                end if
              end if
            end if
          end if
        end do
      end do
    end do
    !
    !--------------------------------------------------------------------
    ! 2. Perform aerosol removal computations
    !  - swith do i,k loop, add chrmbc (the below cloud scavenging rate, s^-1)
    !  - Levin & Schwatz
    !--------------------------------------------------------------------
    !
    if ( ichem == 1 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          mc2mo%rembc(j,i,1) = d_zero
        end do
      end do
      do k = 2 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            mc2mo%rembc(j,i,k) = d_zero
            if ( mc2mo%remrat(j,i,k) > d_zero ) then
              do kk = 1 , k - 1
                qcw = mo2mc%qcn(j,i,k)
                mc2mo%rembc(j,i,k) = mc2mo%rembc(j,i,k) + & ![mm/hr]
                  mc2mo%remrat(j,i,kk) * qcw * &
                              (mo2mc%pfs(j,i,k+1)-mo2mc%pfs(j,i,k))*regrav
              end do
            end if
          end do
        end do
      end do
    end if
    !
    !--------------------------------------------------------------------
    ! 3. Convert the accumlated precipitation to appropriate units for
    !    the surface physics and the output
    !--------------------------------------------------------------------
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        prainx = pptsum(j,i)*dtsec
        if ( prainx > dlowval ) then
          mc2mo%rainnc(j,i) = mc2mo%rainnc(j,i) + prainx
          mc2mo%lsmrnc(j,i) = mc2mo%lsmrnc(j,i) + pptsum(j,i)
          mc2mo%trrate(j,i) = pptsum(j,i)
        end if
      end do
    end do

  contains

#include <wlh.inc>
#include <pfesat.inc>
#include <pfwsat.inc>
#include <clwfromt.inc>

    pure real(rkx) function season_factor(lat) result(sf)
      implicit none
      real(rkx) , intent(in) :: lat
      real(rkx) :: theta , delta
      ! Maximum abs value for the declination angle
      real(rkx) , parameter :: dmax = 0.40910517666747085282_rkx
      ! Different phase in the two emispheres
      if ( lat > d_zero ) then
        theta = twopi*mod(calday+(dayspy*d_half),dayspy)/dayspy
      else
        theta = twopi*calday/dayspy
      end if
      delta = 0.006918_rkx - 0.399912_rkx*cos(theta) + &
              0.070257_rkx*sin(theta) -                &
              0.006758_rkx*cos(2.0_rkx*theta) +        &
              0.000907_rkx*sin(2.0_rkx*theta) -        &
              0.002697_rkx*cos(3.0_rkx*theta) +        &
              0.001480_rkx*sin(3.0_rkx*theta)
      sf = (d_one + delta/dmax)/d_two
    end function season_factor

    subroutine sun_cevap
      implicit none
      integer(ik4) :: i , j
      real(rkx) :: xxlat
      ! cevap minimum seasonal paraneter
      real(rkx) , parameter :: mincevap = 1.0e-5_rkx
      do i = ici1 , ici2
        do j = jci1 , jci2
          xxlat = mo2mc%xlat(j,i)
          xcevap(j,i) = max(cevap(j,i) * (d_one - &
                     (sin(abs(xxlat*90.0_rkx/maxlat)*degrad) * &
                      season_factor(xxlat))), mincevap)
        end do
      end do
    end subroutine sun_cevap

  end subroutine subex

end module mod_micro_subex
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
