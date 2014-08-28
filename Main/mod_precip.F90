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
  real(rk8) , pointer , dimension(:,:) :: pptsum
  real(rk8) , pointer , dimension(:,:) :: psf , rainnc , lsmrnc
  real(rk8) , pointer , dimension(:,:,:,:) :: qx3 , qx2 , qxten
  real(rk8) , pointer , dimension(:,:,:) :: t3 , t2 , tten
  real(rk8) , pointer , dimension(:,:,:) :: p3 , qs3 , rh3 , rho3
  real(rk8) , pointer , dimension(:,:,:) :: radcldf , radlqwc
 
  real(rk8) :: qcth

  real(rk8) , parameter :: uch = d_1000*regrav*secph
  real(rk8) , parameter :: alphaice = d_four

  real(rk8) , public , pointer , dimension(:,:,:) :: fcc , remrat , rembc , totc
  real(rk8) , public , pointer , dimension(:,:) :: qck1 , cgul , rh0 , &
    cevap , caccr
  logical :: lchem = .false.

  public :: allocate_mod_precip , init_precip , pcp , cldfrac , condtq

  contains

  subroutine allocate_mod_precip(ichem)
    implicit none
    integer(ik4) , intent(in) :: ichem
    ! This needs to be saved in SAV file
    call getmem3d(fcc,jci1,jci2,ici1,ici2,1,kz,'pcp:fcc')
    call getmem3d(totc,jci1,jci2,ici1,ici2,1,kz,'pcp:totc')
    ! Those not. Note the external, internal change.
    call getmem2d(qck1,jci1,jci2,ici1,ici2,'pcp:qck1')
    call getmem2d(cgul,jci1,jci2,ici1,ici2,'pcp:cgul')
    call getmem2d(rh0,jci1,jci2,ici1,ici2,'pcp:rh0')
    call getmem2d(cevap,jci1,jci2,ici1,ici2,'pcp:cevap')
    call getmem2d(caccr,jci1,jci2,ici1,ici2,'pcp:caccr')
    call getmem2d(pptsum,jci1,jci2,ici1,ici2,'pcp:pptsum')
    if ( ichem == 1 ) then
      lchem = .true.
      call getmem3d(rembc,jci1,jci2,ici1,ici2,1,kz,'pcp:rembc')
      call getmem3d(remrat,jci1,jci2,ici1,ici2,1,kz,'pcp:remrat')
    end if
  end subroutine allocate_mod_precip

  subroutine init_precip(atmslice,atm,aten,sfs,pptnc,cldfra,cldlwc)
    implicit none
    type(slice) , intent(in) :: atmslice
    type(atmstate) , intent(in) :: atm
    type(atmstate) , intent(in) :: aten
    type(surfstate) , intent(in) :: sfs
    real(rk8) , pointer , dimension(:,:) :: pptnc
    real(rk8) , pointer , dimension(:,:,:) :: cldfra , cldlwc

    call assignpnt(atmslice%tb3d,t3)
    call assignpnt(atmslice%pb3d,p3)
    call assignpnt(atmslice%qxb3d,qx3)
    call assignpnt(atmslice%qsb3d,qs3)
    call assignpnt(atmslice%rhb3d,rh3)
    call assignpnt(atmslice%rhob3d,rho3)
    call assignpnt(atm%t,t2)
    call assignpnt(atm%qx,qx2)
    call assignpnt(aten%t,tten)
    call assignpnt(aten%qx,qxten)
    call assignpnt(sfs%psb,psf)
    call assignpnt(sfs%rainnc,rainnc)
    call assignpnt(pptnc,lsmrnc)
    call assignpnt(cldfra,radcldf)
    call assignpnt(cldlwc,radlqwc)
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
    real(rk8) :: dpovg , afc , ppa , pptacc , pptkm1 , pptmax , &
                pptnew , qcincld , qcleft , qcw , qs , rdevap , &
                rh , rhcs , rho , tcel , thog , tk , prainx
    integer(ik4) :: i , j , k , kk
    !   
    !--------------------------------------------------------------------
    ! 1. Compute the precipitation formed in each layer.
    !    The computations are performed from the top to the surface.
    !    - Auto-conversion: similar to SIMEX (Giorgi and Shields 1999).
    !    - Raindrop Accretion:  Beheng (1994); EQN 6
    !    - Raindrop Evaporation:  Sundqvist (1988); EQN 3.4a
    !--------------------------------------------------------------------
    ! 1a. Perform computations for the top layer (layer 1)
    !   maximum precipation rate (total cloud water/dt)

    thog = d_1000*regrav
    ! precipation accumulated from above
    pptsum(:,:) = d_zero
    if ( lchem ) remrat(:,:,:) = d_zero

    do i = ici1 , ici2
      do j = jci1 , jci2
        afc = fcc(j,i,1)                                     ![frac][avg]
        if ( afc > 0.01D0 ) then !   if there is a cloud
          ! 1aa. Compute temperature and humidities with the adjustments
          !      due to convection.
          tk = t3(j,i,1)                                     ![k][avg]
          tcel = tk - tzero                                  ![C][avg]
          ppa = p3(j,i,1)                                    ![Pa][avg]
          rho = ppa/(rgas*tk)                                ![kg/m3][avg]
          qcw = qx3(j,i,1,iqc)                               ![kg/kg][avg]
          ! 1ab. Calculate the in cloud mixing ratio [kg/kg]
          qcincld = qcw/afc                                  ![kg/kg][cld]
          ! 1ac. Compute the maximum precipation rate
          !      (i.e. total cloud water/dt) [kg/kg/s]
          pptmax = qcw/dt                                    ![kg/kg/s][avg]
          ! 1ad. Implement here the formula for qcth.
          !   - Gultepe & Isaac, J. Clim, 1997, v10 p446 table 4, eq 5
          !   - The factor of 1000 converts from g/kg to kg/kg
          !   - The factor of cgul accounts for the fact that the Gultepe
          !     and Isaac equation is for mean cloud water while qcth is the
          !     theshhold for auto-conversion.
          qcth = cgul(j,i)*(d_10**(-0.489D0+0.0134D0*tcel))*d_r1000 
          ! 1ae. Compute the gridcell average autoconversion [kg/k g/s]
          pptnew = qck1(j,i)*(qcincld-qcth)*afc ! [kg/kg/s][avg]
          pptnew = dmin1(dmax1(pptnew,d_zero),pptmax) ![kg/kg/s][avg]
          if ( pptnew < dlowval ) pptnew = d_zero
          if ( pptnew > d_zero ) then !   New precipitation
            ! 1af. Compute the cloud removal rate (for chemistry) [1/s]
            if (lchem) remrat(j,i,1) = pptnew/qcw
            ! 1ag. Compute the amount of cloud water removed by raindrop
            !      accretion [kg/kg/s].  In the layer where the precipitation
            !      is formed, only half of the precipitation is assumed to
            !      accrete. 1aga. Compute the amount of water remaining in the
            !      cloud [kg/kg]
            qcleft = qcw - pptnew*dt ![kg/kg][avg]
            ! 1agb. Add 1/2 of the new precipitation can accrete.
            pptkm1 = d_half*pptnew/afc*rho*dt ![kg/m3][cld]
            ! 1agc. Accretion [kg/kg/s]=[m3/kg/s]*[kg/kg]*[kg/m3]
            pptacc = caccr(j,i)*qcleft*pptkm1 ![kg/kg/s][avg]
            ! 1agd. Update the precipitation accounting for the
            !       accretion [kg/kg/s]
            pptnew = dmin1(pptmax,pptacc+pptnew) ![kg/kg/s][avg]
            ! 1ah. Accumulate precipitation and convert to kg/m2/s
            dpovg = dsigma(1)*psf(j,i)*thog                  ![kg/m2]
            pptsum(j,i) = pptnew*dpovg ![kg/m2/s][avg]
            ! 1ai. Compute the cloud water tendency [kg/kg/s*cb]
            ! [kg/kg/s*cb][avg]
            qxten(j,i,1,iqc) = qxten(j,i,1,iqc) - pptnew*psf(j,i)
          else  !   Cloud but no new precipitation
            pptsum(j,i) = d_zero ![kg/m2/s][avg]
          end if
        else  !   No cloud
          pptsum(j,i) = d_zero
        end if
      end do
    end do

    ! LAYER TWO TO KZ
    ! 1b. Perform computations for the 2nd layer to the surface
    do k = 2 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          !1ba. Compute temperature and humidities with the adjustments
          !     due to convection.
          tk = t3(j,i,k)                                     ![k][avg]
          tcel = tk - tzero                                  ![C][avg]
          ppa = p3(j,i,k)                                    ![Pa][avg]
          rho = ppa/(rgas*tk)                                ![kg/m3][avg]
          qcw = qx3(j,i,k,iqc)                               ![kg/kg][avg]
          afc = fcc(j,i,k)                                   ![frac][avg]
          qs = pfqsat(tk,ppa)                                ![kg/kg][avg]
          rh = dmin1(dmax1(qx3(j,i,k,iqv)/qs,d_zero),rhmax)  ![frac][avg]
          ! 1bb. Convert accumlated precipitation to kg/kg/s.
          !      Used for raindrop evaporation and accretion.
          dpovg = dsigma(k)*psf(j,i)*thog                    ![kg/m2][avg]
          pptkm1 = pptsum(j,i)/dpovg                         ![kg/kg/s][avg]
          ! 1bc. Compute the raindrop evaporation in the clear portion of
          !      the gridcell.
          !  - It is assumed that raindrops do not evaporate in clouds
          !    and the rainfall from above is evenly distributed in
          !    gridcell (i.e. the gridcell average precipitation is used).
          if ( pptsum(j,i) > d_zero .and. afc < 0.99D0 ) then
            ! 2bca. Compute the clear sky relative humidity
            rhcs = (rh-afc*rhmax)/(d_one-afc)                ![frac][clr]
            rhcs = dmax1(dmin1(rhcs,rhmax),d_zero)           ![frac][clr]
            ! 2bcb. Raindrop evaporation [kg/kg/s]
            rdevap = cevap(j,i)*(rhmax-rhcs)*dsqrt(pptsum(j,i))*(d_one-afc)
                                                           ![kg/kg/s][avg]
            rdevap = dmin1((qs-qx3(j,i,k,iqv))/dt,rdevap)  ![kg/kg/s][avg]
            rdevap = dmin1(dmax1(rdevap,d_zero),pptkm1)    ![kg/kg/s][avg]
            ! 2bcc. Update the precipitation accounting for the raindrop
            !       evaporation [kg/m2/s]
            pptsum(j,i) = pptsum(j,i) - rdevap*dpovg       ![kg/m2/s][avg]
            ! 2bcf. Compute the water vapor tendency [kg/kg/s*cb]
            ![kg/kg/s*cb][avg]
            qxten(j,i,k,iqv) = qxten(j,i,k,iqv) + rdevap*psf(j,i)
            ! 2bcf. Compute the temperature tendency [K/s*cb]
            ![k/s*cb][avg]
            tten(j,i,k) = tten(j,i,k) - wlhvocp*rdevap*psf(j,i)
          else
            !   no precipitation from above
            rdevap = d_zero                                  ![kg/kg/s][avg]
          end if
          ! 1bd. Compute the autoconversion and accretion [kg/kg/s]
          if ( afc > 0.01D0 ) then !   if there is a cloud
            ! 1bda. Calculate the in cloud mixing ratio [kg/kg]
            qcincld = qcw/afc                                ![kg/kg][cld]
            ! 1bdb. Compute the maximum precipation rate
            !       (i.e. total cloud water/dt) [kg/kg/s]
            pptmax = qcw/dt                                  ![kg/kg/s][cld]
            ! 1bdc. Implement the Gultepe & Isaac formula for qcth.
            ![kg/kg][cld]
            qcth = cgul(j,i)*(d_10**(-0.489D0+0.0134D0*tcel))*d_r1000
            ! 1bdd. Compute the gridcell average autoconversion [kg/kg/s]
            pptnew = qck1(j,i)*(qcincld-qcth)*afc            ![kg/kg/s][avg]
            pptnew = dmin1(dmax1(pptnew,d_zero),pptmax)      ![kg/kg/s][avg]
            if ( pptnew < dlowval ) pptnew = d_zero
            ! 1be. Compute the cloud removal rate (for chemistry) [1/s]
            if ( lchem .and. pptnew > d_zero ) remrat(j,i,k) = pptnew/qcw
            ! 1bf. Compute the amount of cloud water removed by raindrop
            !      accretion [kg/kg/s].  In the layer where the precipitation
            !      is formed, only half of the precipitation can accrete.
            if ( pptkm1 > d_zero .or. pptnew > d_zero ) then
              ! 1bfa. Compute the amount of water remaining in the cloud [kg/kg]
              qcleft = qcw-pptnew*dt                         ![kg/kg][avg]
              ! 1bfb. Add 1/2 of the new precipitation to the accumulated
              !       precipitation [kg/m3]
              pptkm1 = (pptkm1+d_half*pptnew/afc)*rho*dt     ![kg/m3][cld]
              ! 1bfc. accretion [kg/kg/s]
              pptacc = caccr(j,i)*qcleft*pptkm1              ![kg/kg/s][avg]
              ! 1bfd. Update the precipitation accounting for the
              !       accretion [kg/kg/s]
              pptnew = dmin1(pptmax,pptacc+pptnew)           ![kg/kg/s][avg]
            end if
            ! 1bg. Accumulate precipitation and convert to kg/m2/s
            pptsum(j,i) = pptsum(j,i) + pptnew*dpovg        ![kg/m2/s][avg]
            ! 1bh. Compute the cloud water tendency [kg/kg/s*cb]
            ![kg/kg/s*cb][avg]
            qxten(j,i,k,iqc) = qxten(j,i,k,iqc) - pptnew*psf(j,i)
          else
            pptnew = d_zero ![kg/kg/s][avg]
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
          rembc(j,i,1) = d_zero
        end do
      end do
      do k = 2 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            rembc(j,i,k) = d_zero
            if ( remrat(j,i,k) > d_zero ) then
              do kk = 1 , k - 1
                rembc(j,i,k) = rembc(j,i,k) + remrat(j,i,kk)*qx3(j,i,kk,iqc) * &
                             psf(j,i)*dsigma(kk)*uch          ![mm/hr]
              end do
              ! the below cloud precipitation rate is now used
              ! directly in chemistry
!             rembc(j,i,k) = 6.5D0*1.0D-5*rembc(j,i,k)**0.68D0   ![s^-1]
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
    real(rk8) :: exlwc , rh0adj , tcel
    integer(ik4) :: i , j , k
    !--- add hanzhenyu 20140520
    real(rk8) :: botm , rm , excld , qcld

    exlwc = d_zero
    botm = d_zero
    rm = d_zero
    excld = d_zero
    qcld = d_zero
    !------

    !-----------------------------------------
    ! 1.  Determine large-scale cloud fraction
    !-----------------------------------------
    do k = 1 , kz
      ! Adjusted relative humidity threshold
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( t3(j,i,k) >= tc0 ) then
            rh0adj = rh0(j,i)
          else ! high cloud (less subgrid variability)
            rh0adj = rhmax-(rhmax-rh0(j,i))/(d_one+0.15D0*(tc0-t3(j,i,k)))
          end if
!         if (ipptls == 1) then
            if ( rh3(j,i,k) >= rhmax) then        ! full cloud cover
              fcc(j,i,k) = d_one
            else if ( rh3(j,i,k) <= rh0adj ) then  ! no cloud cover
              fcc(j,i,k) = d_zero
            else                                         ! partial cloud cover
              fcc(j,i,k) = d_one-dsqrt(d_one-(rh3(j,i,k)-rh0adj) / &
                           (rhmax-rh0adj))
              ! if (ipptls == 1) totc(j,i,k) = qx3(j,i,k,iqc)
              ! if (ipptls > 1)  totc(j,i,k) = qx3(j,i,k,iqc) + qx3(j,i,k,iqi)
              fcc(j,i,k) = dmin1(dmax1(fcc(j,i,k),0.01D0),0.99D0)
             end if
!         else
!           if ( rh3(j,i,k) >= d_one) then        ! full cloud cover
!             fcc(j,i,k) = d_one
!           else 
!             fcc(j,i,k) = (rh3(j,i,k)**0.25D0)* &
!                   (d_one-dexp((-100.0D0*(qx3(j,i,k,iqc)+qx3(j,i,k,iqi))/ & 
!                   ((d_one-rh3(j,i,k))*qs3(j,i,k))**0.49D0)))
!             fcc(j,i,k) = dmin1(dmax1(fcc(j,i,k),0.01D0),0.99D0)
!           end if !  rh0 threshold
!           Test CF either 1 or 0
!           if (qx3(j,i,k,iqc)+qx3(j,i,k,iqi)>minqx) then
!             fcc(j,i,k)=d_one
!           else 
!             fcc(j,i,k)=d_zero
!           end if 
!         end if
          !----------------------------------------------------------------
          ! Correction:
          !   Ivan Guettler, 14.10.2010.
          ! Based on: Vavrus, S. and Waliser D., 2008, 
          ! An Improved Parameterization for Simulating Arctic Cloud Amount
          ! in the CCSM3 Climate Model, J. Climate 
          !----------------------------------------------------------------
          if ( p3(j,i,k) >= 75000.0D0 ) then
            ! Clouds below 750hPa
            if ( qx3(j,i,k,iqv) <= 0.003D0 ) then
              fcc(j,i,k) = fcc(j,i,k) * &
                     dmax1(0.15D0,dmin1(d_one,qx3(j,i,k,iqv)/0.003D0))
              !
              ! Tuğba Öztürk mod for Siberia
              !
              ! fcc(j,i,k) = (rh3(j,i,k)**0.25D0)* &
              !      (d_one-dexp((-100.0D0*qx3(j,i,k,iqc)) / &
              !     ((d_one-rh3(j,i,k))*qs3(j,i,k))**0.49D0))
              !
              !
            end if
          end if
          !----------------------------------------------------------------
          ! End of the correction.
          !----------------------------------------------------------------
        end do
      end do
    end do
    !-----------------------------------------------------------------
    ! 2.  Combine large-scale and convective fraction and liquid water
    !     to be passed into radiation.
    !-----------------------------------------------------------------
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          ! Calculate total condensate, in ipptls = 1 it is given by the
          ! liquid only, in ipptls = 2 it is given by liquid and ice
          if ( ipptls == 1 ) then
            totc(j,i,k) = qx3(j,i,k,iqc)
          else if ( ipptls == 2 ) then
            totc(j,i,k) = qx3(j,i,k,iqc) + alphaice*qx3(j,i,k,iqi)
          end if
          !-------- add hanzhenyu 20140520
          if ( iconvlwp == 2 ) then
            ! To calculate cloud fraction using the semi-empirical formula
            ! of Xu and Randall (1996, JAS)
            qcld = totc(j,i,k)
            if ( qcld < 1.0D-12 ) then    ! no cloud cover
              excld = d_zero
            else if ( rh3(j,i,k) >= d_one ) then  ! full cloud cover
              excld = d_one
            else
              botm = exp( 0.49D0*log((d_one-rh3(j,i,k))*qs3(j,i,k)) )
              rm = exp(0.25D0*log(rh3(j,i,k)))
              if ( 100.D0*(qcld/botm) > 25.0D0 ) then
                excld = min(0.99D0, max(0.01D0 , rm))
              else
                excld = min(0.99D0, max(0.01D0 , &
                        rm*(d_one-exp(-100.D0*(qcld/botm)))))  
              end if
              ! xu and randall
            end if
            !------ for low cloud
            ! a) clouds below 750hPa, extremely cold conditions,
            !     when no cld microphy
            if ( ipptls == 1 .and. p3(j,i,k) >= 75000.0D0 .and. &
                    qx3(j,i,k,iqv) <= 0.003D0 ) then
              excld = fcc(j,i,k)   ! reduced rh-cld
            end if
            !--------
            ! Cloud Water Volume    
            ! kg gq / kg dry air * kg dry air / m3 * 1000 = g qc /  m3
            if ( qcld > minqx ) then
              exlwc = qcld*rho3(j,i,k)*d_1000
            else
              exlwc = 0.001D0
            end if
            !------------
          else if ( iconvlwp == 1 .or. iconvlwp == 0 ) then
            !------------
            ! Adjusted relative humidity threshold
            if ( t3(j,i,k) >= tc0 ) then
              rh0adj = rh0(j,i)
            else ! high cloud (less subgrid variability)
              rh0adj = rhmax-(rhmax-rh0(j,i))/(d_one+0.15D0*(tc0-t3(j,i,k)))
            end if
            if ( rh3(j,i,k) >= rhmax ) then        ! full cloud cover
              excld = d_one
            else if ( rh3(j,i,k) <= rh0adj ) then  ! no cloud cover
              excld = d_zero
            else                                   ! partial cloud cover
              excld = d_one-dsqrt(d_one-(rh3(j,i,k)-rh0adj) / &
                            (rhmax-rh0adj))
              excld = dmin1(dmax1(excld,0.01D0),0.99D0)
            end if !  rh0 threshold
            !----------------------------------------------------------------
            ! Based on: Vavrus, S. and Waliser D., 2008, 
            !----------------------------------------------------------------
            ! Clouds below 750hPa
            if ( p3(j,i,k) >= 75000.0D0 .and. qx3(j,i,k,iqv) <= 0.003D0 ) then
              excld = excld * dmax1(0.15D0,dmin1(d_one,qx3(j,i,k,iqv)/0.003D0))
            end if
            !------------------- end hanzhenyu 
            !--------
            if ( totc(j,i,k) > minqx ) then
              if ( iconvlwp == 1 ) then
                ! Apply the parameterisation based on temperature to the
                ! the large scale clouds.
                tcel = t3(j,i,k)-tzero
                if ( tcel < -50D0 ) then
                  exlwc = 0.001D0
                else
                  exlwc = 0.127D+00 + 6.78D-03 * tcel + &
                          1.29D-04 * tcel**2 + 8.36D-07 * tcel**3
                  if ( exlwc > 0.300D0) exlwc = 0.300D0
                  if ( exlwc < 0.001D0) exlwc = 0.001D0
                end if
              else
                ! Cloud Water Volume
                ! kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
                exlwc = totc(j,i,k)*rho3(j,i,k)*d_1000
              end if
            else
              exlwc = 0.001D0
            end if
          end if   ! end if of iconvlwp , hanzy
          radlqwc(j,i,k) = (radcldf(j,i,k)*radlqwc(j,i,k) + excld*exlwc) / &
               dmax1(radcldf(j,i,k)+excld,0.01D0)
          radcldf(j,i,k) = dmin1(dmax1(radcldf(j,i,k),excld),cftotmax)
          if ( iconvlwp == 2 .and. ipptls == 2 ) then
            fcc(j,i,k) = excld
          end if
!hanzy    radlqwc(j,i,k) = (radcldf(j,i,k)*radlqwc(j,i,k) + &
!hanzy    fcc(j,i,k)*exlwc) / max(radcldf(j,i,k)+fcc(j,i,k),0.01D0)
!hanzy    radcldf(j,i,k) = min(max(fcc(j,i,k),radcldf(j,i,k)),cftotmax)
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
  subroutine condtq(psc)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(in) :: psc
    !
    ! rhc    - Relative humidity at ktau+1
    ! rh0adj - Adjusted relative humidity threshold at ktau+1
    ! fccc   - Cloud fraction at ktau+1
    !
    real(rk8) :: qccs , qvcs , tmp1 , tmp2 , tmp3
    real(rk8) :: dqv , exces , fccc , pres , qvc_cld , qvs , &
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
            write(stderr,*) 'At global J : ',global_dot_jstart+j
            write(stderr,*) 'At global I : ',global_dot_istart+i
            write(stderr,*) 'At global K : ',k
          end if
#endif
          qvcs = dmax1((qx2(j,i,k,iqv)+dt*qxten(j,i,k,iqv))/psc(j,i),minqx)
          qccs = dmax1((qx2(j,i,k,iqc)+dt*qxten(j,i,k,iqc))/psc(j,i),d_zero)
          !-----------------------------------------------------------
          !     2.  Compute the cloud condensation/evaporation term.
          !-----------------------------------------------------------
          ! 2a. Calculate the saturation mixing ratio and relative humidity
          pres = (psc(j,i)*hsigma(k)+ptop)*1000.0D0
          qvs = pfqsat(tmp3,pres)
          rhc = dmax1(qvcs/qvs,dlowval)

          r1 = d_one/(d_one+wlhv*wlhv*qvs/(rwat*cpd*tmp3*tmp3))
     
          ! 2b. Compute the relative humidity threshold at ktau+1
          if ( tmp3 > tc0 ) then
            rh0adj = rh0(j,i)
          else ! high cloud (less subgrid variability)
            rh0adj = rhmax - (rhmax-rh0(j,i))/(d_one+0.15D0*(tc0-tmp3))
          end if
     
          ! 2c. Compute the water vapor in excess of saturation
          if ( rhc >= rhmax .or. rhc < rh0adj ) then ! Full or no cloud cover
            dqv = qvcs - qvs*conf ! Water vapor in excess of sat
            tmp1 = r1*dqv
            fccc=d_zero
          else                                       ! Partial cloud cover
            fccc = d_one - dsqrt(d_one-(rhc-rh0adj)/(rhmax-rh0adj))
            !---- add hanzhenyu 20140520
            if ( pres >= 75000.0D0 .and. qvcs <= 0.003D0 ) then
              ! Clouds below 750hPa
              fccc = fccc * dmax1(0.15D0,dmin1(d_one,qvcs/0.003D0))
            end if
            !--------------------
            fccc = dmin1(dmax1(fccc,0.01D0),d_one)
            qvc_cld = dmax1((qs3(j,i,k)+dt*qxten(j,i,k,iqv)/psc(j,i)),d_zero)
            dqv = qvc_cld - qvs*conf  ! qv diff between predicted qv_c
            tmp1 = r1*dqv*fccc        ! grid cell average
          end if
     
          ! 2d. Compute the new cloud water + old cloud water
          exces = qccs + tmp1
          if ( exces >= d_zero ) then ! Some cloud is left
            tmp2 = tmp1/dt
          else                        ! The cloud evaporates
            tmp2 = -qccs/dt
          end if
          !-----------------------------------------------------------
          !     3.  Compute the tendencies.
          !-----------------------------------------------------------
          qxten(j,i,k,iqv) = qxten(j,i,k,iqv) - psc(j,i)*tmp2
          qxten(j,i,k,iqc) = qxten(j,i,k,iqc) + psc(j,i)*tmp2
          tten(j,i,k) = tten(j,i,k) + psc(j,i)*tmp2*wlhvocp
        end do
      end do
    end do
  end subroutine condtq

end module mod_precip
