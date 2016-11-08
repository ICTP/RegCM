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

module mod_precip
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
  use mod_humid , only : clwfromt

  implicit none

  private
  !
  ! Precip sum beginning from top
  !
  real(rkx) , pointer , dimension(:,:) :: pptsum
  integer(ik4) , pointer , dimension(:,:) :: ldmsk
  real(rkx) , pointer , dimension(:,:) :: psb , psc , rainnc , lsmrnc
  real(rkx) , pointer , dimension(:,:) :: th700 , sfps
  real(rkx) , pointer , dimension(:,:,:) :: xqc , xqi
  real(rkx) , pointer , dimension(:,:,:,:) :: qx3 , qx2 , qxten , chi3
  real(rkx) , pointer , dimension(:,:,:) :: t3 , t2 , tten , th3
  real(rkx) , pointer , dimension(:,:,:) :: p3 , p2 , qs3 , rh3 , rho3
  real(rkx) , pointer , dimension(:,:,:) :: radcldf , radlqwc
  real(rkx) , pointer , dimension(:,:,:) :: pfcc
  real(rkx) , pointer , dimension(:,:,:) :: premrat
  real(rkx) , pointer , dimension(:,:,:) :: prembc
  real(rkx) , pointer , dimension(:,:,:) :: totc , pccn
  real(rkx) , pointer , dimension(:,:,:) :: dia_qcr , dia_qcl , dia_acr

  real(rkx) :: maxlat

  real(rkx) , parameter :: thog = d_1000*regrav
  real(rkx) , parameter :: uch = thog*secph
  real(rkx) , parameter :: alphaice = d_four

  real(rkx) , public , pointer , dimension(:,:) :: qck1 , cgul , rh0 , &
    cevap , xcevap , caccr
  real(rkx) , public , pointer , dimension(:,:,:) :: dqc
  real(rkx) , public , pointer , dimension(:,:,:) :: qvn , qln

  logical :: l_lat_hack = .false.
  public :: allocate_mod_precip , init_precip , pcp , cldfrac , condtq

  real(rkx) , parameter :: rhow = 1000.0_rkx
  real(rkx) , parameter :: qvmin = 1.0e-8_rkx
  real(rkx) , parameter :: pptmin = 1.0e-20_rkx

  contains

#include <pfesat.inc>
#include <pfwsat.inc>

  subroutine allocate_mod_precip
    implicit none
    ! Those not. Note the external, internal change.
    call getmem2d(qck1,jci1,jci2,ici1,ici2,'pcp:qck1')
    call getmem2d(cgul,jci1,jci2,ici1,ici2,'pcp:cgul')
    call getmem2d(rh0,jci1,jci2,ici1,ici2,'pcp:rh0')
    call getmem2d(cevap,jci1,jci2,ici1,ici2,'pcp:cevap')
    if ( l_lat_hack ) then
      call getmem2d(xcevap,jci1,jci2,ici1,ici2,'pcp:xcevap')
    else
      call assignpnt(cevap,xcevap)
    end if
    call getmem2d(caccr,jci1,jci2,ici1,ici2,'pcp:caccr')
    call getmem2d(pptsum,jci1,jci2,ici1,ici2,'pcp:pptsum')
    call getmem3d(dqc,jci1,jci2,ici1,ici2,1,kz,'pcp:dqc')
    call getmem3d(totc,jci1,jci2,ici1,ici2,1,kz,'pcp:totc')
    if ( ipptls == 2 ) then
      call getmem3d(qln,jci1,jci2,ici1,ici2,1,kz,'pcp:qln')
    end if
  end subroutine allocate_mod_precip

  subroutine init_precip
    use mod_atm_interface , only : mddom , atms , atm2 , aten , sfs , &
                                   pptnc , cldfra , cldlwc , fcc ,    &
                                   remrat , rembc , qdiag , ccn
    use mod_mppparam , only : maxall
    implicit none
    call maxall(maxval(mddom%xlat),maxlat)
    call assignpnt(mddom%ldmsk,ldmsk)
    call assignpnt(atms%tb3d,t3)
    call assignpnt(atms%th3d,th3)
    call assignpnt(atms%th700,th700)
    call assignpnt(atms%pb3d,p3)
    call assignpnt(atms%qxb3d,qx3)
    call assignpnt(atms%qxb3d,qvn,iqv)
    call assignpnt(atms%qxb3d,xqc,iqc)
    call assignpnt(atms%qsb3d,qs3)
    call assignpnt(atms%rhb3d,rh3)
    call assignpnt(atms%rhob3d,rho3)
    call assignpnt(atms%ps2d,sfps)
    call assignpnt(atm2%t,t2)
    call assignpnt(atm2%qx,qx2)
    call assignpnt(atm2%pr,p2)
    call assignpnt(aten%t,tten)
    call assignpnt(aten%qx,qxten)
    call assignpnt(sfs%psb,psb)
    call assignpnt(sfs%psc,psc)
    call assignpnt(sfs%rainnc,rainnc)
    call assignpnt(pptnc,lsmrnc)
    call assignpnt(cldfra,radcldf)
    call assignpnt(cldlwc,radlqwc)
    call assignpnt(fcc,pfcc)
    if ( ichem == 1 ) then
      call assignpnt(atms%chib3d,chi3)
      call assignpnt(remrat,premrat)
      call assignpnt(rembc,prembc)
      call assignpnt(ccn,pccn)
      if ( idiag == 1 ) then
        call assignpnt(qdiag%qcr,dia_qcr)
        call assignpnt(qdiag%qcl,dia_qcl)
        call assignpnt(qdiag%acr,dia_acr)
      end if
    end if
    if ( ipptls == 1 ) then
      call assignpnt(atms%qxb3d,qln,iqc)
    else if ( ipptls == 2 ) then
      call assignpnt(atms%qxb3d,xqi,iqi)
    end if
  end subroutine init_precip
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
  subroutine pcp
    implicit none
    real(rkx) :: dpovg , afc , pptacc , pptkm1 , pptmax ,       &
                pptnew , qcleft , qcw , qs , rdevap , qcincl ,  &
                rh , rhcs , rho , tcel , prainx , qcth
    integer(ik4) :: i , j , k , kk
    logical :: lsecind
    !
    ! 0. Compute dqc
    !
    lsecind = (ichem == 1 .and. iaerosol == 1 .and. iindirect == 2)
    if ( l_lat_hack ) then
      call sun_cevap
    end if
    if ( ipptls == 2 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            qln(j,i,k) = xqc(j,i,k) + xqi(j,i,k)
          end do
        end do
      end do
    end if
    if ( lsecind ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! include aerosol second indirect effect on threshold
            ! auto-conversion
            ! rcrit is a critical cloud radius for cloud
            ! water undergoing autoconversion
            ! pccn = number of ccn /m3
            ! In cloud mixing ratio [kg/kg]
            qcincl = qln(j,i,k)/pfcc(j,i,k)
            qcth = pccn(j,i,k)*(4.0_rkx/3.0_rkx)*mathpi * &
                          ((rcrit*1e-6_rkx)**3)*rhow
            if ( idiag == 1 ) then
              dia_qcr(j,i,k) = qcth
              dia_qcl(j,i,k) = qcincl
            end if
            dqc(j,i,k) = qcincl - qcth
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
            tcel = t3(j,i,k) - tzero   ![C][avg]
            ! In cloud mixing ratio [kg/kg]
            qcincl = qln(j,i,k)/pfcc(j,i,k)
            qcth = cgul(j,i)*(d_10**(-0.489_rkx+0.0134_rkx*tcel))*d_r1000
            dqc(j,i,k) = qcincl - qcth
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
      premrat(:,:,:) = d_zero
      if ( lsecind .and. idiag == 1 ) then
        dia_acr(:,:,:) = d_zero
      end if
    end if
    ! 1a. Perform computations for the top layer (layer 1)
    !   maximum precipation rate (total cloud water/dt)
    do i = ici1 , ici2
      do j = jci1 , jci2
        afc = pfcc(j,i,1)   ![frac][avg]
        qcw = qln(j,i,1)    ![kg/kg][avg]
        pptnew = d_zero
        if ( afc > lowcld ) then ! if there is a cloud
          ! 1ac. Compute the maximum precipation rate
          !      (i.e. total cloud water/dt) [kg/kg/s]
          pptmax = max(qcw,d_zero)/dt             ![kg/kg/s][avg]
          ! 1ae. Compute the gridcell average autoconversion [kg/k g/s]
          pptnew = qck1(j,i)*dqc(j,i,1)*afc       ![kg/kg/s][avg]
          pptnew = min(max(pptnew,d_zero),pptmax) ![kg/kg/s][avg]
          if ( pptnew > pptmin ) then !   New precipitation
            ! 1ag. Compute the amount of cloud water removed by raindrop
            !      accretion [kg/kg/s].  In the layer where the precipitation
            !      is formed, only half of the precipitation is assumed to
            !      accrete. 1aga. Compute the amount of water remaining in the
            !      cloud [kg/kg]
            qcleft = qcw - pptnew*dt ![kg/kg][avg]
            ! 1agb. Add 1/2 of the new precipitation can accrete.
            rho = rho3(j,i,1)                  ![kg/m3][avg]
            pptkm1 = d_half*pptnew/afc*rho*dt  ![kg/m3][cld]
            ! 1agc. Accretion [kg/kg/s]=[m3/kg/s]*[kg/kg]*[kg/m3]
            pptacc = caccr(j,i)*qcleft*pptkm1  ![kg/kg/s][avg]
            ! 1agd. Update the precipitation accounting for the
            !       accretion [kg/kg/s]
            pptnew = min(pptmax,pptacc+pptnew) ![kg/kg/s][avg]
          end if
        end if
        ! 1ah. Accumulate precipitation and convert to kg/m2/s
        if ( pptnew > pptmin ) then
          dpovg = dsigma(1)*psb(j,i)*thog    ![kg/m2]
          pptsum(j,i) = pptnew*dpovg         ![kg/m2/s][avg]
          ! 1ai. Compute the cloud water tendency [kg/kg/s*cb]
          ! [kg/kg/s*cb][avg]
          qxten(j,i,1,iqc) = qxten(j,i,1,iqc) - pptnew*psb(j,i)
          ! 1af. Compute the cloud removal rate (for chemistry) [1/s]
          if ( ichem == 1 ) then
            premrat(j,i,1) = pptnew/qcw
            if ( lsecind .and. idiag == 1 ) then
              dia_acr(j,i,1) = pptnew
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
          pptnew = d_zero
          ! 1bb. Convert accumlated precipitation to kg/kg/s.
          !      Used for raindrop evaporation and accretion.
          dpovg = dsigma(k)*psb(j,i)*thog                    ![kg/m2][avg]
          if ( pptsum(j,i) > d_zero ) then
            pptkm1 = pptsum(j,i)/dpovg                       ![kg/kg/s][avg]
          else
            pptkm1 = d_zero
          end if
          ! 1bc. Compute the raindrop evaporation in the clear portion of
          !      the gridcell.
          !  - It is assumed that raindrops do not evaporate in clouds
          !    and the rainfall from above is evenly distributed in
          !    gridcell (i.e. the gridcell average precipitation is used).
          afc = pfcc(j,i,k)                              ![frac][avg]
          if ( pptkm1 > pptmin .and. afc < hicld ) then
            ! 2bca. Compute the clear sky relative humidity
            rh = rh3(j,i,k)                              ![frac][avg]
            rhcs = (rh-afc*rhmax)/(hicld-afc)            ![frac][clr]
            rhcs = max(min(rhcs,rhmax),rhmin)            ![frac][clr]
            ! 2bcb. Raindrop evaporation [kg/kg/s]
            rdevap = xcevap(j,i)*(rhmax-rhcs)*sqrt(pptsum(j,i))*(hicld-afc)
            qs = pfwsat(t3(j,i,k),p3(j,i,k))             ![kg/kg][avg]
            rdevap = min((qs-qvn(j,i,k))/dt,rdevap)      ![kg/kg/s][avg]
            rdevap = min(max(rdevap,d_zero),pptkm1)      ![kg/kg/s][avg]
            ! 2bcc. Update the precipitation accounting for the raindrop
            !       evaporation [kg/m2/s]
            pptsum(j,i) = max(pptsum(j,i)-rdevap*dpovg,d_zero) ![kg/m2/s][avg]
            ! 2bcf. Compute the water vapor tendency [kg/kg/s*cb]
            ![kg/kg/s*cb][avg]
            qxten(j,i,k,iqv) = qxten(j,i,k,iqv) + rdevap*psb(j,i)
            ! 2bcf. Compute the temperature tendency [K/s*cb]
            ![k/s*cb][avg]
            tten(j,i,k) = tten(j,i,k) - wlhvocp*rdevap*psb(j,i)
          end if
          qcw = qln(j,i,k)   ![kg/kg][avg]
          ! 1bd. Compute the autoconversion and accretion [kg/kg/s]
          if ( afc > lowcld ) then ! if there is a cloud
            ! 1bdb. Compute the maximum precipation rate
            !       (i.e. total cloud water/dt) [kg/kg/s]
            pptmax = max(qcw,d_zero)/dt                  ![kg/kg/s][avg]
            ! 1bdd. Compute the gridcell average autoconversion [kg/kg/s]
            pptnew = qck1(j,i)*dqc(j,i,k)*afc            ![kg/kg/s][avg]
            pptnew = min(max(pptnew,d_zero),pptmax)      ![kg/kg/s][avg]
            ! 1bf. Compute the amount of cloud water removed by raindrop
            !      accretion [kg/kg/s].  In the layer where the precipitation
            !      is formed, only half of the precipitation can accrete.
            if ( pptkm1 > pptmin .or. pptnew > pptmin ) then
              ! 1bfa. Compute the amount of water remaining in the cloud [kg/kg]
              qcleft = qcw-pptnew*dt                         ![kg/kg][avg]
              ! 1bfb. Add 1/2 of the new precipitation to the accumulated
              !       precipitation [kg/m3]
              rho = rho3(j,i,k)                              ![kg/m3][avg]
              pptkm1 = (pptkm1+d_half*pptnew/afc)*rho*dt     ![kg/m3][cld]
              ! 1bfc. accretion [kg/kg/s]
              pptacc = caccr(j,i)*qcleft*pptkm1              ![kg/kg/s][avg]
              ! 1bfd. Update the precipitation accounting for the
              !       accretion [kg/kg/s]
              pptnew = min(max(pptacc+pptnew,d_zero),pptmax)
            end if
          end if
          if ( pptnew > pptmin ) then
            ! 1bg. Accumulate precipitation and convert to kg/m2/s
            pptsum(j,i) = pptsum(j,i) + pptnew*dpovg       ![kg/m2/s][avg]
            ! 1bh. Compute the cloud water tendency [kg/kg/s*cb]
            ![kg/kg/s*cb][avg]
            qxten(j,i,k,iqc) = qxten(j,i,k,iqc) - pptnew*psb(j,i)
            ! 1be. Compute the cloud removal rate (for chemistry) [1/s]
            if ( ichem == 1 ) then
              premrat(j,i,k) = pptnew/qcw
              if (lsecind .and. idiag == 1 ) then
                dia_acr(j,i,k) = pptnew
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
          prembc(j,i,1) = d_zero
        end do
      end do
      do k = 2 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            prembc(j,i,k) = d_zero
            if ( premrat(j,i,k) > d_zero ) then
              do kk = 1 , k - 1
                qcw = qln(j,i,k)
                prembc(j,i,k) = prembc(j,i,k) + & ![mm/hr]
                  premrat(j,i,kk) * qcw * psb(j,i) * dsigma(kk) * uch
              end do
              ! the below cloud precipitation rate is now used
              ! directly in chemistry
!             prembc(j,i,k) = 6.5_rkx*1.0e-5_rkx * |
!                        prembc(j,i,k)**0.68_rkx   ![s^-1]
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
          rainnc(j,i) = rainnc(j,i) + prainx
          lsmrnc(j,i) = lsmrnc(j,i) + pptsum(j,i)
        end if
      end do
    end do

  contains

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
              0.070257_rkx*sin(theta) -              &
              0.006758_rkx*cos(2.0_rkx*theta) +        &
              0.000907_rkx*sin(2.0_rkx*theta) -        &
              0.002697_rkx*cos(3.0_rkx*theta) +        &
              0.001480_rkx*sin(3.0_rkx*theta)
      sf = (d_one + delta/dmax)/d_two
    end function season_factor

    subroutine sun_cevap
      use mod_atm_interface , only : mddom
      implicit none
      integer(ik4) :: i , j
      real(rkx) :: xxlat
      ! cevap minimum seasonal paraneter
      real(rkx) , parameter :: mincevap = 1.0e-5_rkx
      do i = ici1 , ici2
        do j = jci1 , jci2
          xxlat = mddom%xlat(j,i)
          xcevap(j,i) = max(cevap(j,i) * (d_one - &
                     (sin(abs(xxlat*90.0_rkx/maxlat)*degrad) * &
                      season_factor(xxlat))), mincevap)
        end do
      end do
    end subroutine sun_cevap

  end subroutine pcp
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                 c
  ! This subroutine computes the fractional cloud coverage and      c
  ! liquid water content (in cloud value).  Both are use in         c
  ! radiation.                                                      c
  !                                                                 c
  ! The fractional coverage of large scale clouds is a function of  c
  ! relative humidity, using the relationship of sundqvist et       c
  ! al., 1989.  The relative humidity at which clouds begin to      c
  ! form is lower over land than ocean, due to the greater number   c
  ! of cloud condensation nucleii.                                  c
  !                                                                 c
  ! The fracional coverage of convective clouds is passed in from   c
  ! the convection scheme.                                          c
  !                                                                 c
  ! The large-scale and convective clouds are combined as follows:  c
  ! 1) If the convective cloud fraction > large scale fraction, the c
  ! convective fraction and water content are used (this occurs     c
  ! infrequently).                                                  c
  ! 2) Otherwise, the cloud fraction equals the large-scale         c
  ! fraction AND the water content is a weighted average of both    c
  ! types.                                                          c
  !                                                                 c
  ! Note: the incloud water content (g/m3) is passed to radiation   c
  !                                                                 c
  ! See Pal et al (2000) for more info.                             c
  !                                                                 c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine cldfrac
    implicit none
    real(rkx) :: exlwc , rh0adj
    integer(ik4) :: i , j , k
    real(rkx) :: pres , botm , rm , qcld

    if ( ipptls == 2 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            totc(j,i,k) = (xqc(j,i,k) + alphaice*xqi(j,i,k)) / &
                                (d_one+qvn(j,i,k))
          end do
        end do
      end do
    else
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            totc(j,i,k) = xqc(j,i,k)/(d_one+qvn(j,i,k))
          end do
        end do
      end do
    end if

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------

    if ( icldfrac == 1 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! To calculate cloud fraction using the semi-empirical formula
            ! of Xu and Randall (1996, JAS)
            if ( rh3(j,i,k) >= rhmax ) then  ! full cloud cover
              pfcc(j,i,k) = hicld
            else if ( rh3(j,i,k) <= rhmin ) then
              pfcc(j,i,k) = lowcld
            else
              qcld = totc(j,i,k)
              botm = exp( 0.49_rkx*log((rhmax-rh3(j,i,k))*qs3(j,i,k)) )
              rm = exp(0.25_rkx*log(rh3(j,i,k)))
              if ( 100._rkx*(qcld/botm) > 25.0_rkx ) then
                pfcc(j,i,k) = min(hicld, max(lowcld , rm))
              else
                pfcc(j,i,k) = min(hicld, max(lowcld , &
                        rm*(d_one-exp(-100._rkx*(qcld/botm)))))
              end if
            end if
          end do
        end do
      end do
    else
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! Adjusted relative humidity threshold
            pres = p3(j,i,k)
            if ( t3(j,i,k) > tc0 ) then
              rh0adj = rh0(j,i)
            else ! high cloud (less subgrid variability)
              rh0adj = rhmax-(rhmax-rh0(j,i))/(d_one+0.15_rkx*(tc0-t3(j,i,k)))
            end if
            rh0adj = max(rhmin,min(rh0adj,rhmax))
            if ( rh3(j,i,k) >= rhmax ) then     ! full cloud cover
              pfcc(j,i,k) = hicld
            else if ( rh3(j,i,k) <= rhmin ) then
              pfcc(j,i,k) = lowcld
            else
              ! Use Sundqvist (1989) formula
              pfcc(j,i,k) = d_one-sqrt((rhmax-rh3(j,i,k)) / &
                                       (rhmax-rh0adj))
              pfcc(j,i,k) = min(max(pfcc(j,i,k),lowcld),hicld)
            end if
          end do
        end do
      end do
    end if
    !
    ! Correction:
    !   Ivan Guettler, 14.10.2010.
    ! Based on: Vavrus, S. and Waliser D., 2008,
    ! An Improved Parameterization for Simulating Arctic Cloud Amount
    ! in the CCSM3 Climate Model, J. Climate
    !
    if ( ipptls == 1 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! clouds below 750hPa, extremely cold conditions,
            !  when no cld microphy
            if ( p3(j,i,k) >= 75000.0_rkx .and. &
                 qvn(j,i,k) <= 0.003_rkx ) then
              pfcc(j,i,k) = pfcc(j,i,k) * &
                     max(0.15_rkx,min(d_one,qvn(j,i,k)/0.003_rkx))
              !
              ! Tuğba Öztürk mod for Siberia
              !
              ! pfcc(j,i,k) = (rh3(j,i,k)**0.25_rkx)* &
              !      (d_one-dexp((-100.0_rkx*qx3(j,i,k,iqc)) / &
              !     ((rhmax-rh3(j,i,k))*qs3(j,i,k))**0.49_rkx))
              !
              !
            end if
          end do
        end do
      end do
    end if

    !------------------------------------------
    ! 1a. Determine Marine stratocumulus clouds
    !------------------------------------------

    if ( icldmstrat == 1 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( ldmsk(j,i) == 0 ) then
              if ( p3(j,i,k) >= 70000.0_rkx ) then
                ! Klein, S. A., and D. L. Hartmann,
                ! The seasonal cycle of low stratiform clouds,
                ! J. Climate, 6, 1587-1606, 1993
                pfcc(j,i,k) = min(hicld,max(pfcc(j,i,k), &
                      (th700(j,i)-th3(j,i,k)) * 0.057_rkx - 0.5573_rkx))
              end if
            end if
          end do
        end do
      end do
    end if

    !-----------------------------------------------------------------
    ! 2.  Combine large-scale and convective fraction and liquid water
    !     to be passed into radiation.
    !-----------------------------------------------------------------

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          exlwc = d_zero
          ! Cloud Water Volume
          ! kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
          if ( iconvlwp == 1 ) then
            ! Apply the parameterisation based on temperature to the
            ! the large scale clouds.
            if ( pfcc(j,i,k) > lowcld ) then
              exlwc = clwfromt(t3(j,i,k))
            end if
          else
            ! NOTE : IN CLOUD HERE IS NEEDED !!!
            if ( pfcc(j,i,k) > lowcld ) then
              ! In g / m^3
              exlwc = ((totc(j,i,k)*d_1000)/pfcc(j,i,k))*rho3(j,i,k)
            end if
          end if
          if ( radcldf(j,i,k) > lowcld .or. pfcc(j,i,k) > lowcld ) then
            ! get maximum cloud fraction between cumulus and large scale
            radlqwc(j,i,k) = (exlwc * pfcc(j,i,k) + &
                              radlqwc(j,i,k) * radcldf(j,i,k)) / &
                              (radcldf(j,i,k) + pfcc(j,i,k))
            radcldf(j,i,k) = max(radcldf(j,i,k),pfcc(j,i,k))
          else
            radcldf(j,i,k) = lowcld
            radlqwc(j,i,k) = d_zero
          end if
          radcldf(j,i,k) = min(max(radcldf(j,i,k),lowcld),cftotmax)
        end do
      end do
    end do
  end subroutine cldfrac
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                 c
  ! This subroutine computes the condensational or evaporational    c
  ! heating term and the fallout term of precipitation from the     c
  ! explicit moisture scheme.                                       c
  !                                                                 c
  ! ---the condensational or evaporational term are one step        c
  !    adjustment based on asai (1965, j. meteo. soc. japan).       c
  !                                                                 c
  ! ---modified to include the effects of partial cloud cover       c
  !    (see Pal et al 2000).  When partial clouds exist, the qxten  c
  !    in/out of the clear and cloudy portions of the grid cell is  c
  !    assumed to be at the same rate (i.e., if there is 80% cloud  c
  !    cover, .2 of qxten goes to raising qv in the clear region    c
  !    and .8 goes to condensation or evaporation of qc in the      c
  !    cloudy portion).                                             c
  !                                                                 c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine condtq
    implicit none
    !
    ! rhc    - Relative humidity at ktau+1
    ! rh0adj - Adjusted relative humidity threshold at ktau+1
    ! fccc   - Cloud fraction at ktau+1
    !
    real(rkx) :: qccs , qvcs , tmp1 , tmp2 , tmp3
    real(rkx) :: dqv , exces , fccc , pres , qvc_cld , qvs , &
               r1 , rh0adj , rhc
    integer(ik4) :: i , j , k

    !---------------------------------------------------------------------
    !     1.  Compute t, qv, and qc at tau+1 without condensational term
    !---------------------------------------------------------------------
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          tmp3 = (t2(j,i,k)+dt*tten(j,i,k))/psc(j,i)
#ifdef DEBUG
          if ( tmp3 < d_zero ) then
            write(stderr,*) 'Time ktau = ', ktau
            write(stderr,*) 'Consistency TEMPERATURE ERROR in condtq (T < 0K)'
            write(stderr,*) 'At global J : ',j
            write(stderr,*) 'At global I : ',i
            write(stderr,*) 'At global K : ',k
          end if
#endif
          qvcs = max((qx2(j,i,k,iqv)+dt*qxten(j,i,k,iqv)),qvmin)/psc(j,i)
          qccs = max((qx2(j,i,k,iqc)+dt*qxten(j,i,k,iqc)),d_zero)/psc(j,i)
          !-----------------------------------------------------------
          !     2.  Compute the cloud condensation/evaporation term.
          !-----------------------------------------------------------
          ! 2a. Calculate the saturation mixing ratio and relative humidity
          pres = p2(j,i,k)
          qvs = pfwsat(tmp3,pres)
          rhc = min(max(qvcs/qvs,rhmin),rhmax)

          r1 = d_one/(d_one+wlhv*wlhv*qvs/(rwat*cpd*tmp3*tmp3))

          ! 2b. Compute the relative humidity threshold at ktau+1
          if ( tmp3 > tc0 ) then
            rh0adj = rh0(j,i)
          else ! high cloud (less subgrid variability)
            rh0adj = rhmax - (rhmax-rh0(j,i))/(d_one+0.15_rkx*(tc0-tmp3))
          end if
          rh0adj = max(rhmin,min(rh0adj,rhmax))
          if ( rhc < rh0adj ) then      ! No cloud cover
            dqv = qvcs - qvs*conf
          else if ( rhc >= rhmax ) then ! Full or no cloud cover
            dqv = qvcs - qvs*conf
          else
            fccc = d_one - sqrt((rhmax-rhc)/(rhmax-rh0adj))
            if ( pres >= 75000.0_rkx .and. qvcs <= 0.003_rkx ) then
              fccc = fccc * max(0.15_rkx, min(d_one,qvcs/0.003_rkx))
            end if
            fccc = min(max(fccc,d_zero),d_one)
            qvc_cld = max((qs3(j,i,k)+dt*qxten(j,i,k,iqv)/psc(j,i)),d_zero)
            dqv = fccc * (qvc_cld - qvs*conf)  ! qv diff between predicted qv_c
          end if

          ! 2c. Compute the water vapor in excess of saturation
          tmp1 = r1*dqv               ! grid cell average

          ! 2d. Compute the new cloud water + old cloud water
          exces = qccs + tmp1
          if ( exces >= d_zero ) then ! Some cloud is left
            tmp2 = tmp1/dt
          else                        ! The cloud evaporates
            tmp2 = -(qccs*conf)/dt
          end if
          !-----------------------------------------------------------
          !     3.  Compute the tendencies.
          !-----------------------------------------------------------
          if ( abs(tmp2) > dlowval ) then
            qxten(j,i,k,iqv) = qxten(j,i,k,iqv) - psc(j,i)*tmp2
            qxten(j,i,k,iqc) = qxten(j,i,k,iqc) + psc(j,i)*tmp2
            tten(j,i,k) = tten(j,i,k) + psc(j,i)*tmp2*wlhvocp
          end if
        end do
      end do
    end do
  end subroutine condtq

end module mod_precip
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
