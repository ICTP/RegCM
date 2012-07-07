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
  use mod_runparams
  use mod_memutil
  use mod_atm_interface , only : atmstate , slice , surfstate
!
  private
!
! Precip sum beginning from top
  real(dp) , pointer , dimension(:,:) :: pptsum
  real(dp) , pointer , dimension(:,:) :: psf , rainnc , lsmrnc
  real(dp) , pointer , dimension(:,:,:,:) :: qx3 , qx2 , qxten
  real(dp) , pointer , dimension(:,:,:) :: t3 , t2 , tten
  real(dp) , pointer , dimension(:,:,:) :: p3 , qs3 , rh3 , rho3
  real(dp) , pointer , dimension(:,:,:) :: cldfra , cldlwc
 
  real(dp) :: qcth , aprdiv
!
  real(dp) , parameter :: uch = d_1000*regrav*secph
  real(dp) , parameter :: lowq = 1.0D-8
!
  real(dp) , public , pointer , dimension(:,:,:) :: fcc , remrat , rembc
  real(dp) , public , pointer , dimension(:,:) :: qck1 , cgul , rh0
  real(dp) , public :: caccr , cevap , rhmax , tc0 , fcmax , conf
  ! TAO 2/8/11:
  ! Flag for using convective liquid water path as the large-scale
  ! liquid water path (iconvlwp=1)
  integer , public :: iconvlwp
  logical :: lchem = .false.
!
  public :: allocate_mod_precip , init_precip , pcp , cldfrac , condtq
!
  contains
!
    subroutine allocate_mod_precip(ichem)
      implicit none
      integer , intent(in) :: ichem
      ! This needs to be saved in SAV file
      call getmem3d(fcc,jce1,jce2,ice1,ice2,1,kz,'pcp:fcc')
      ! Those not. Note the external, internal change.
      call getmem2d(qck1,jci1,jci2,ici1,ici2,'pcp:qck1')
      call getmem2d(cgul,jci1,jci2,ici1,ici2,'pcp:cgul')
      call getmem2d(rh0,jci1,jci2,ici1,ici2,'pcp:rh0')
      if ( ichem == 1 ) then
        lchem = .true.
        call getmem3d(rembc,jci1,jci2,ici1,ici2,1,kz,'pcp:rembc')
        call getmem3d(remrat,jci1,jci2,ici1,ici2,1,kz,'pcp:remrat')
      end if
    end subroutine allocate_mod_precip
!
    subroutine init_precip(atmslice,atm,aten,sfs,pptnc,radcldf,radlqwc)
      implicit none
      type(slice) , intent(in) :: atmslice
      type(atmstate) , intent(in) :: atm
      type(atmstate) , intent(in) :: aten
      type(surfstate) , intent(in) :: sfs
      real(dp) , pointer , dimension(:,:) :: pptnc
      real(dp) , pointer , dimension(:,:,:) :: radcldf , radlqwc

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
      call assignpnt(radcldf,cldfra)
      call assignpnt(radlqwc,cldlwc)

      call getmem2d(pptsum,jci1,jci2,ici1,ici2,'pcp:pptsum')

      aprdiv = d_one/dble(ntsrf)
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
!
    real(dp) :: dpovg , es , afc , ppa , pptacc , pptkm1 , pptmax ,   &
               pptnew , qcincld , qcleft , qcw , qs , rdevap , &
               rh , rhcs , rho , tcel , thog , tk , prainx
    integer :: i , j , k , kk
!
!   
!--------------------------------------------------------------------
!   1. Compute the precipitation formed in each layer.
!      The computations are performed from the top to the surface.
!      - Auto-conversion: similar to SIMEX (Giorgi and Shields 1999).
!      - Raindrop Accretion:  Beheng (1994); EQN 6
!      - Raindrop Evaporation:  Sundqvist (1988); EQN 3.4a
!--------------------------------------------------------------------
   
!   1a. Perform computations for the top layer (layer 1)
!   maximum precipation rate (total cloud water/dt)

    thog = d_1000*regrav
!   precipation accumulated from above
    pptsum(:,:) = d_zero
    if ( lchem ) remrat(:,:,:) = d_zero

    do i = ici1 , ici2
      do j = jci1 , jci2
        afc = fcc(j,i,1)                                     ![frac][avg]
        if ( afc > 0.01D0 ) then !   if there is a cloud
!       1aa. Compute temperature and humidities with the adjustments
!            due to convection.
          tk = t3(j,i,1)                                     ![k][avg]
          tcel = tk - tzero                                  ![C][avg]
          ppa = p3(j,i,1)*d_1000                             ![Pa][avg]
          rho = ppa/(rgas*tk)                                ![kg/m3][avg]
          qcw = qx3(j,i,1,iqc)                               ![kg/kg][avg]
!         1ab. Calculate the in cloud mixing ratio [kg/kg]
          qcincld = qcw/afc                                  ![kg/kg][cld]
!         1ac. Compute the maximum precipation rate
!              (i.e. total cloud water/dt) [kg/kg/s]
          pptmax = qcw/dt                                    ![kg/kg/s][avg]
!         1ad. Implement here the formula for qcth.
!              - Gultepe & Isaac, J. Clim, 1997, v10 p446 table 4, eq 5
!              - The factor of 1000 converts from g/kg to kg/kg
!              - The factor of cgul accounts for the fact that the Gultepe
!                and Isaac equation is for mean cloud water while qcth is the
!                theshhold for auto-conversion.
          qcth = cgul(j,i)*(d_10**(-0.489D0+0.0134D0*tcel))*d_r1000 
                                                             ![kg/kg][cld]
!         1ae. Compute the gridcell average autoconversion [kg/k g/s]
          pptnew = qck1(j,i)*(qcincld-qcth)*afc              ![kg/kg/s][avg]
          pptnew = dmin1(dmax1(pptnew,d_zero),pptmax)        ![kg/kg/s][avg]
          if ( pptnew < dlowval ) pptnew = d_zero

          if ( pptnew > d_zero ) then !   New precipitation
!           1af. Compute the cloud removal rate (for chemistry) [1/s]
            if (lchem) remrat(j,i,1) = pptnew/qcw
!           1ag. Compute the amount of cloud water removed by raindrop
!                accretion [kg/kg/s].  In the layer where the precipitation
!                is formed, only half of the precipitation is assumed to
!                accrete. 1aga. Compute the amount of water remaining in the
!                cloud [kg/kg]
            qcleft = qcw - pptnew*dt                         ![kg/kg][avg]
!           1agb. Add 1/2 of the new precipitation can accrete.
            pptkm1 = d_half*pptnew/afc*rho*dt                ![kg/m3][cld]
!           1agc. Accretion [kg/kg/s]=[m3/kg/s]*[kg/kg]*[kg/m3]
            pptacc = caccr*qcleft*pptkm1                     ![kg/kg/s][avg]
!           1agd. Update the precipitation accounting for the accretion
!           [kg/kg/s]
            pptnew = dmin1(pptmax,pptacc+pptnew)             ![kg/kg/s][avg]
!           1ah. Accumulate precipitation and convert to kg/m2/s
            dpovg = dsigma(1)*psf(j,i)*thog                  ![kg/m2]
            pptsum(j,i) = pptnew*dpovg                         ![kg/m2/s][avg]
!           1ai. Compute the cloud water tendency [kg/kg/s*cb]
            qxten(j,i,1,iqc) = qxten(j,i,1,iqc) - pptnew*psf(j,i) ![kg/kg/s*cb][avg]
          else  !   Cloud but no new precipitation
            pptsum(j,i) = d_zero                               ![kg/m2/s][avg]
          end if
        else  !   No cloud
          pptsum(j,i) = d_zero                                 ![kg/m2/s][avg]
        end if
      end do
    end do

!     LAYER TWO TO KZ
!     1b. Perform computations for the 2nd layer to the surface
    do k = 2 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2

!         1ba. Compute temperature and humidities with the adjustments
!              due to convection.
          tk = t3(j,i,k)                                     ![k][avg]
          tcel = tk - tzero                                  ![C][avg]
          ppa = p3(j,i,k)*d_1000                             ![Pa][avg]
          rho = ppa/(rgas*tk)                                ![kg/m3][avg]
          qcw = qx3(j,i,k,iqc)                               ![kg/kg][avg]
          afc = fcc(j,i,k)                                   ![frac][avg]
          if ( tcel > d_zero ) then
            es = svp1*d_1000*dexp(svp2*tcel/(tk-svp3))       ![Pa][avg]
          else
            es = svp4*d_1000*dexp(svp5-svp6/tk)              ![Pa][avg]
          end if
          qs = ep2*es/(ppa-es)                               ![kg/kg][avg]
          rh = dmin1(dmax1(qx3(j,i,k,iqv)/qs,d_zero),rhmax)  ![frac][avg]
    
!         1bb. Convert accumlated precipitation to kg/kg/s.
!              Used for raindrop evaporation and accretion.
          dpovg = dsigma(k)*psf(j,i)*thog                    ![kg/m2][avg]
          pptkm1 = pptsum(j,i)/dpovg                         ![kg/kg/s][avg]
   
!         1bc. Compute the raindrop evaporation in the clear portion of
!              the gridcell.
!            - It is assumed that raindrops do not evaporate in clouds
!              and the rainfall from above is evenly distributed in
!                gridcell (i.e. the gridcell average precipitation is used).
          if ( pptsum(j,i) > d_zero .and. afc < 0.99D0 ) then
!           2bca. Compute the clear sky relative humidity
            rhcs = (rh-afc*rhmax)/(d_one-afc)                ![frac][clr]
            rhcs = dmax1(dmin1(rhcs,rhmax),d_zero)           ![frac][clr]
!           2bcb. Raindrop evaporation [kg/kg/s]
            rdevap = cevap*(rhmax-rhcs)*dsqrt(pptsum(j,i))*(d_one-afc)
                                                           ![kg/kg/s][avg]
            rdevap = dmin1((qs-qx3(j,i,k,iqv))/dt,rdevap)  ![kg/kg/s][avg]
            rdevap = dmin1(dmax1(rdevap,d_zero),pptkm1)    ![kg/kg/s][avg]
!           2bcc. Update the precipitation accounting for the raindrop
!                 evaporation [kg/m2/s]
            pptsum(j,i) = pptsum(j,i) - rdevap*dpovg       ![kg/m2/s][avg]
!           2bcf. Compute the water vapor tendency [kg/kg/s*cb]
            qxten(j,i,k,iqv) = qxten(j,i,k,iqv) + rdevap*psf(j,i)
                                                           ![kg/kg/s*cb][avg]
!           2bcf. Compute the temperature tendency [K/s*cb]
            tten(j,i,k) = tten(j,i,k) - wlhvocp*rdevap*psf(j,i)
                                                           ![k/s*cb][avg]
          else
            !   no precipitation from above
            rdevap = d_zero                                  ![kg/kg/s][avg]
          end if
    
!         1bd. Compute the autoconversion and accretion [kg/kg/s]
          if ( afc > 0.01D0 ) then !   if there is a cloud
!           1bda. Calculate the in cloud mixing ratio [kg/kg]
            qcincld = qcw/afc                                ![kg/kg][cld]
!           1bdb. Compute the maximum precipation rate
!                 (i.e. total cloud water/dt) [kg/kg/s]
            pptmax = qcw/dt                                  ![kg/kg/s][cld]
!           1bdc. Implement the Gultepe & Isaac formula for qcth.
            qcth = cgul(j,i)*(d_10**(-0.489D0+0.0134D0*tcel))*d_r1000
                                                           ![kg/kg][cld]
!           1bdd. Compute the gridcell average autoconversion [kg/kg/s]
            pptnew = qck1(j,i)*(qcincld-qcth)*afc            ![kg/kg/s][avg]
            pptnew = dmin1(dmax1(pptnew,d_zero),pptmax)      ![kg/kg/s][avg]
            if ( pptnew < dlowval ) pptnew = d_zero
!           1be. Compute the cloud removal rate (for chemistry) [1/s]
            if ( lchem .and. pptnew > d_zero ) remrat(j,i,k) = pptnew/qcw
    
!           1bf. Compute the amount of cloud water removed by raindrop
!                accretion [kg/kg/s].  In the layer where the precipitation
!                is formed, only half of the precipitation can accrete.
            if ( pptkm1 > d_zero .or. pptnew > d_zero ) then
!             1bfa. Compute the amount of water remaining in the cloud [kg/kg]
              qcleft = qcw-pptnew*dt                         ![kg/kg][avg]
!             1bfb. Add 1/2 of the new precipitation to the accumulated
!                   precipitation [kg/m3]
              pptkm1 = (pptkm1+d_half*pptnew/afc)*rho*dt     ![kg/m3][cld]
!             1bfc. accretion [kg/kg/s]
              pptacc = caccr*qcleft*pptkm1                   ![kg/kg/s][avg]
!             1bfd. Update the precipitation accounting for the
!                   accretion [kg/kg/s]
              pptnew = dmin1(pptmax,pptacc+pptnew)           ![kg/kg/s][avg]
            end if
!           1bg. Accumulate precipitation and convert to kg/m2/s
            pptsum(j,i) = pptsum(j,i) + pptnew*dpovg        ![kg/m2/s][avg]
!           1bh. Compute the cloud water tendency [kg/kg/s*cb]
            qxten(j,i,k,iqc) = qxten(j,i,k,iqc) - pptnew*psf(j,i) ![kg/kg/s*cb][avg]
          else
            pptnew = d_zero                                  ![kg/kg/s][avg]
          end if
        end do
      end do
    end do
!     
!--------------------------------------------------------------------
!     2. Perform aerosol removal computations
!       - swith do i,k loop, add chrmbc (the below cloud scavenging rate, s^-1)
!       - Levin & Schwatz
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
!   the below cloud precipitation rate is now used directly in chemistry
!             rembc(j,i,k) = 6.5D0*1.0D-5*rembc(j,i,k)**0.68D0   ![s^-1]

            end if
          end do
        end do
      end do
    end if
!     
!--------------------------------------------------------------------
!     3. Convert the accumlated precipitation to appropriate units for
!        the surface physics and the output
!--------------------------------------------------------------------
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        prainx = pptsum(j,i)*dtsec
        if ( prainx > dlowval ) then
          rainnc(j,i) = rainnc(j,i) + prainx
          if ( ktau == 0 .and. debug_level > 2 ) then
            lsmrnc(j,i) = lsmrnc(j,i) + pptsum(j,i)
          else
            lsmrnc(j,i) = lsmrnc(j,i) + pptsum(j,i)*aprdiv
          end if
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
!
    implicit none
!
    real(dp) :: exlwc , rh0adj
    integer :: i , j , k
!
!--------------------------------------------------------------------
!     1.  Determine large-scale cloud fraction
!--------------------------------------------------------------------
    do k = 1 , kz
      ! Adjusted relative humidity threshold
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( t3(j,i,k) > tc0 ) then
            rh0adj = rh0(j,i)
          else ! high cloud (less subgrid variability)
            rh0adj = rhmax-(rhmax-rh0(j,i))/(d_one+0.15D0*(tc0-t3(j,i,k)))
          end if
          if ( rh3(j,i,k) >= rhmax ) then        ! full cloud cover
            fcc(j,i,k) = d_one
          else if ( rh3(j,i,k) <= rh0adj ) then  ! no cloud cover
            fcc(j,i,k) = d_zero
          else                                   ! partial cloud cover
            fcc(j,i,k) = d_one-dsqrt(d_one-(rh3(j,i,k)-rh0adj) / &
                          (rhmax-rh0adj))
            fcc(j,i,k) = dmin1(dmax1(fcc(j,i,k),0.01D0),0.99D0)
          end if !  rh0 threshold
!---------------------------------------------------------------------
!         Correction:
!         Ivan Guettler, 14.10.2010.
!         Based on: Vavrus, S. and Waliser D., 2008, 
!         An Improved Parameterization for Simulating Arctic Cloud Amount
!         in the CCSM3 Climate Model, J. Climate 
!---------------------------------------------------------------------
          if ( p3(j,i,k) >= 75.0D0 ) then
            ! Clouds below 750hPa
            if ( qx3(j,i,k,iqv) <= 0.003D0 ) then
              fcc(j,i,k) = fcc(j,i,k) * &
                     dmax1(0.15D0,dmin1(d_one,qx3(j,i,k,iqv)/0.003D0))
            end if
          end if
!---------------------------------------------------------------------
!         End of the correction.
!---------------------------------------------------------------------
        end do
      end do
    end do

!--------------------------------------------------------------------
!     2.  Combine large-scale and convective fraction and liquid water
!         to be passed into radiation.
!--------------------------------------------------------------------
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          ! Cloud Water Volume
          ! kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
          exlwc = qx3(j,i,k,iqc)*rho3(j,i,k)*d_1000

          ! temperature dependance for convective cloud water content
          ! in g/m3 (Lemus et al., 1997)
          cldlwc(j,i,k)  = 0.127D+00 + 6.78D-03*(t3(j,i,k)-tzero) + &
                         1.29D-04* (t3(j,i,k)-tzero)**d_two  +    &
                         8.36D-07*(t3(j,i,k)-tzero)**d_three

          if ( cldlwc(j,i,k) > 0.3D+00 ) cldlwc(j,i,k) = 0.3D+00
          if ( (t3(j,i,k)-tzero) < -50D+00 ) cldlwc(j,i,k) = 0.001D+00
          ! Apply the parameterisation based on temperature to the
          ! convective fraction AND the large scale clouds :
          ! the large scale cloud water content is not really used by
          ! current radiation code, needs further evaluation.
          !TAO: but only apply this parameterization to large scale LWC 
          !if the user specifies it
          if (iconvlwp == 1) exlwc = cldlwc(j,i,k)
          cldlwc(j,i,k) = (cldfra(j,i,k)*cldlwc(j,i,k)+fcc(j,i,k)*exlwc) / &
                        dmax1(cldfra(j,i,k)+fcc(j,i,k),0.01D0)
          cldfra(j,i,k) = dmin1(dmax1(cldfra(j,i,k),fcc(j,i,k)),fcmax)
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
!
    implicit none
!
    real(dp) , pointer , dimension(:,:) , intent(in) :: psc
!
!   rhc    - Relative humidity at ktau+1
!   rh0adj - Adjusted relative humidity threshold at ktau+1
!   fccc   - Cloud fraction at ktau+1
!
    real(dp) :: qccs , qvcs , tmp1 , tmp2 , tmp3
    real(dp) :: dqv , exces , fccc , pres , qvc_cld , qvs , &
               r1 , rh0adj , rhc , satvp
    integer :: i , j , k

    !---------------------------------------------------------------------
    !     1.  Compute t, qv, and qc at tau+1 without condensational term
    !---------------------------------------------------------------------
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          tmp3 = (t2(j,i,k)+dt*tten(j,i,k))/psc(j,i)
          qvcs = dmax1((qx2(j,i,k,iqv)+dt*qxten(j,i,k,iqv))/psc(j,i),lowq)
          qccs = dmax1((qx2(j,i,k,iqc)+dt*qxten(j,i,k,iqc))/psc(j,i),dlowval)
          !-----------------------------------------------------------
          !     2.  Compute the cloud condensation/evaporation term.
          !-----------------------------------------------------------
          ! 2a. Calculate the saturation mixing ratio and relative humidity
          pres = (hsigma(k)*psc(j,i)+ptop)*d_1000
          if ( tmp3 > tzero ) then
            satvp = svp1*d_1000*dexp(svp2*(tmp3-tzero)/(tmp3-svp3))
          else
            satvp = svp4*d_1000*dexp(svp5-svp6/tmp3)
          end if
          qvs = dmax1(ep2*satvp/(pres-satvp),lowq)
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
          else                                       ! Partial cloud cover
            fccc = d_one - dsqrt(d_one-(rhc-rh0adj)/(rhmax-rh0adj))
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
!
end module mod_precip
