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
  use mod_atm_interface , only : atmstate , slice , surfpstate , surfstate
!
  private
!
! Precip sum beginning from top
  real(8) , pointer , dimension(:) :: pptsum
  real(8) , pointer , dimension(:,:) :: psf , rainnc , lsmrnc
  real(8) , pointer , dimension(:,:,:) :: t3 , p3 , qv3 , qc3 , qs3 , rh3 , rho3
  real(8) , pointer , dimension(:,:,:) :: t2 , qc2 , qv2
  real(8) , pointer , dimension(:,:,:) :: tten , qvten , qcten
  real(8) , pointer , dimension(:,:) :: cldfra , cldlwc
 

  real(8) :: qcth
  integer :: istart , istopx
  real(8) , pointer , dimension(:,:) :: qccs , tmp1 , tmp2 , tmp3
!
  real(8) , parameter :: uch = d_1000*regrav*secph
  real(8) , parameter :: lowq = 1.0D-30
!
  real(8) , public , pointer , dimension(:,:,:) :: fcc
  real(8) , public , pointer , dimension(:,:) :: qck1 , cgul , rh0,remrat,rembc
  real(8) , public :: caccr , cevap , rhmax , tc0 , fcmax , conf
  ! TAO 2/8/11:
  ! Flag for using convective liquid water path as the large-scale
  ! liquid water path (iconvlwp=1)
  integer , public :: iconvlwp
!
  public :: allocate_mod_precip , init_precip , pcp , cldfrac , condtq
!
  contains
!
    subroutine allocate_mod_precip
      implicit none
      call getmem3d(fcc,1,iy,1,kz,1,jxp,'pcp:fcc')
      call getmem2d(qck1,1,iy,1,jxp,'pcp:qck1')
      call getmem2d(cgul,1,iy,1,jxp,'pcp:cgul')
      call getmem2d(rh0,1,iy,1,jxp,'pcp:rh0')
      call getmem2d(remrat,1,iy,1,kz,'pcp:remrat')
      call getmem2d(rembc,1,iy,1,kz,'pcp:rembc')
    end subroutine allocate_mod_precip
!
    subroutine init_precip(atmslice,atm,tendency,sps,surface,lsmrainnc, &
                           radcldf,radlqwc)
      implicit none
      type(slice) , intent(in) :: atmslice
      type(atmstate) , intent(in) :: atm
      type(atmstate) , intent(in) :: tendency
      type(surfpstate) , intent(in) :: sps
      type(surfstate) , intent(in) :: surface
      real(8) , pointer , dimension(:,:) :: lsmrainnc
      real(8) , pointer , dimension(:,:) :: remrat , rembc
      real(8) , pointer , dimension(:,:) :: radcldf , radlqwc

      call assignpnt(atmslice%tb3d,t3)
      call assignpnt(atmslice%pb3d,p3)
      call assignpnt(atmslice%qvb3d,qv3)
      call assignpnt(atmslice%qcb3d,qc3)
      call assignpnt(atmslice%qsb3d,qs3)
      call assignpnt(atmslice%rhb3d,rh3)
      call assignpnt(atmslice%rhob3d,rho3)
      call assignpnt(atm%t,t2)
      call assignpnt(atm%qv,qv2)
      call assignpnt(atm%qc,qc2)
      call assignpnt(tendency%t,tten)
      call assignpnt(tendency%qc,qcten)
      call assignpnt(tendency%qv,qvten)
      call assignpnt(sps%ps,psf)
      call assignpnt(surface%rainnc,rainnc)
      call assignpnt(lsmrainnc,lsmrnc)
 !     call assignpnt(remrat,chrmrat)
 !     call assignpnt(rembc,chrmbc)
      call assignpnt(radcldf,cldfra)
      call assignpnt(radlqwc,cldlwc)

      istart = lbound(rainnc,1) + 1
      istopx = ubound(rainnc,1) - 2

      call getmem1d(pptsum,istart,istopx,'pcp:pptsum')
      call getmem2d(qccs,istart,istopx,1,kz,'pcp:qccs')
      call getmem2d(tmp1,istart,istopx,1,kz,'pcp:tmp1')
      call getmem2d(tmp2,istart,istopx,1,kz,'pcp:tmp2')
      call getmem2d(tmp3,istart,istopx,1,kz,'pcp:tmp3')
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
    subroutine pcp(jstart,jstop)

    implicit none
!
    integer , intent(in) :: jstart , jstop
!
    real(8) :: aprdiv , dpovg , es , afc , p , pptacc ,               &
               pptkm1 , pptmax , pptnew , q , qcincld , qcleft , qcw ,&
               qs , rdevap , rh , rhcs , rho , tcel , thog , tk ,     &
               uconv , prainx
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
    pptsum(:) = d_zero
    remrat(istart:istopx,1:kz) = d_zero

    do j = jstart , jstop
!
      do i = istart , istopx
     
        afc = fcc(i,1,j)                                     ![frac][avg]
     
        if ( afc > 0.01D0 ) then !   if there is a cloud
!       1aa. Compute temperature and humidities with the adjustments
!            due to convection.
          q = qv3(i,1,j)                                     ![kg/kg][avg]
          tk = t3(i,1,j)                                     ![k][avg]
          tcel = tk - tzero                                  ![C][avg]
          p = p3(i,1,j)*d_1000                               ![Pa][avg]
          rho = p/(rgas*tk)                                  ![kg/m3][avg]
          qcw = qc3(i,1,j)                                   ![kg/kg][avg]
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
          qcth = cgul(i,j)*(d_10**(-0.489D0+0.0134D0*tcel))*d_r1000 
                                                             ![kg/kg][cld]
!         1ae. Compute the gridcell average autoconversion [kg/k g/s]
          pptnew = qck1(i,j)*(qcincld-qcth)*afc              ![kg/kg/s][avg]
          pptnew = dmin1(dmax1(pptnew,d_zero),pptmax)        ![kg/kg/s][avg]
          if ( pptnew < dlowval ) pptnew = d_zero

          if ( pptnew > d_zero ) then !   New precipitation
!           1af. Compute the cloud removal rate (for chemistry) [1/s]
            remrat(i,1) = pptnew/qcw
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
            dpovg = dsigma(1)*psf(i,j)*thog                  ![kg/m2]
            pptsum(i) = pptnew*dpovg                         ![kg/m2/s][avg]
!           1ai. Compute the cloud water tendency [kg/kg/s*cb]
            qcten(i,1,j) = qcten(i,1,j) - pptnew*psf(i,j)    ![kg/kg/s*cb][avg]
          else  !   Cloud but no new precipitation
            pptsum(i) = d_zero                               ![kg/m2/s][avg]
          end if
        else  !   No cloud
          pptsum(i) = d_zero                                 ![kg/m2/s][avg]
        end if
     
      end do
     
!     LAYER TWO TO KZ
!     1b. Perform computations for the 2nd layer to the surface
      do k = 2 , kz
        do i = istart , istopx
     
!         1ba. Compute temperature and humidities with the adjustments
!              due to convection.
          q = qv3(i,k,j)                                     ![kg/kg][avg]
          tk = t3(i,k,j)                                     ![k][avg]
          tcel = tk - tzero                                  ![C][avg]
          p = p3(i,k,j)*d_1000                               ![Pa][avg]
          rho = p/(rgas*tk)                                  ![kg/m3][avg]
          qcw = qc3(i,k,j)                                   ![kg/kg][avg]
          afc = fcc(i,k,j)                                   ![frac][avg]
          if ( tcel > d_zero ) then
            es = svp1*d_1000*dexp(svp2*tcel/(tk-svp3))       ![Pa][avg]
          else
            es = svp4*d_1000*dexp(svp5-svp6/tk)              ![Pa][avg]
          end if
          qs = ep2*es/(p-es)                                 ![kg/kg][avg]
          rh = dmin1(dmax1(q/qs,d_zero),rhmax)               ![frac][avg]
     
!         1bb. Convert accumlated precipitation to kg/kg/s.
!              Used for raindrop evaporation and accretion.
          dpovg = dsigma(k)*psf(i,j)*thog                    ![kg/m2][avg]
          pptkm1 = pptsum(i)/dpovg                           ![kg/kg/s][avg]
     
!         1bc. Compute the raindrop evaporation in the clear portion of
!              the gridcell.
!              - It is assumed that raindrops do not evaporate in clouds
!                and the rainfall from above is evenly distributed in
!                gridcell (i.e. the gridcell average precipitation is used).
          if ( pptsum(i) > d_zero .and. afc < 0.99D0 ) then
!           2bca. Compute the clear sky relative humidity
            rhcs = (rh-afc*rhmax)/(d_one-afc)                ![frac][clr]
            rhcs = dmax1(dmin1(rhcs,rhmax),d_zero)           ![frac][clr]
!           2bcb. Raindrop evaporation [kg/kg/s]
            rdevap = cevap*(rhmax-rhcs)*dsqrt(pptsum(i))*(d_one-afc)
                                                             ![kg/kg/s][avg]
            rdevap = dmin1((qs-q)/dt,rdevap)                 ![kg/kg/s][avg]
            rdevap = dmin1(dmax1(rdevap,d_zero),pptkm1)      ![kg/kg/s][avg]
!           2bcc. Update the precipitation accounting for the raindrop
!                 evaporation [kg/m2/s]
            pptsum(i) = pptsum(i) - rdevap*dpovg             ![kg/m2/s][avg]
!           2bcf. Compute the water vapor tendency [kg/kg/s*cb]
            qvten(i,k,j) = qvten(i,k,j) + rdevap*psf(i,j)    ![kg/kg/s*cb][avg]
!           2bcf. Compute the temperature tendency [K/s*cb]
            tten(i,k,j) = tten(i,k,j) - wlhvocp*rdevap*psf(i,j)
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
            qcth = cgul(i,j)*(d_10**(-0.489D0+0.0134D0*tcel))*d_r1000
                                                             ![kg/kg][cld]
!           1bdd. Compute the gridcell average autoconversion [kg/kg/s]
            pptnew = qck1(i,j)*(qcincld-qcth)*afc            ![kg/kg/s][avg]
            pptnew = dmin1(dmax1(pptnew,d_zero),pptmax)      ![kg/kg/s][avg]
            if ( pptnew < dlowval ) pptnew = d_zero
!           1be. Compute the cloud removal rate (for chemistry) [1/s]
            if ( pptnew > d_zero ) remrat(i,k) = pptnew/qcw
     
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
            pptsum(i) = pptsum(i) + pptnew*dpovg             ![kg/m2/s][avg]
!           1bh. Compute the cloud water tendency [kg/kg/s*cb]
            qcten(i,k,j) = qcten(i,k,j) - pptnew*psf(i,j)    ![kg/kg/s*cb][avg]
          else
            pptnew = d_zero                                  ![kg/kg/s][avg]
          end if
     
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
        do i = istart , istopx
          rembc(i,1) = d_zero
          do k = 2 , kz
            rembc(i,k) = d_zero
            if ( remrat(i,k) > d_zero ) then
              do kk = 1 , k - 1
                rembc(i,k) = rembc(i,k) + remrat(i,kk)*qc3(i,kk,j) *  &
                             psf(i,j)*dsigma(kk)*uch          ![mm/hr]
              end do
! the below cloud precipitation rate is now used directly in chemistry
!              rembc(i,k) = 6.5D0*1.0D-5*rembc(i,k)**0.68D0   ![s^-1]

            end if
          end do
        end do
      end if
!     
!--------------------------------------------------------------------
!     3. Convert the accumlated precipitation to appropriate units for
!        the surface physics and the output
!--------------------------------------------------------------------
!
      uconv = dtsec
      aprdiv = d_one/dble(ntsrf)
      if ( ktau == 0 ) aprdiv = uconv
      do i = istart , istopx
        prainx = pptsum(i)*uconv
        if ( prainx > dlowval ) then
          rainnc(i,j) = rainnc(i,j) + prainx
          lsmrnc(i,j) = lsmrnc(i,j) + pptsum(i)*aprdiv
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
    subroutine cldfrac(jstart,jstop)
!
    implicit none
!
    integer , intent(in) :: jstart , jstop
!
    real(8) :: exlwc , rh0adj
    integer :: i , j , k
!
    do j = jstart , jstop
!--------------------------------------------------------------------
!     1.  Determine large-scale cloud fraction
!--------------------------------------------------------------------
      do k = 1 , kz
        ! Adjusted relative humidity threshold
        do i = 2 , iym2
          if ( t3(i,k,j) > tc0 ) then
            rh0adj = rh0(i,j)
          else ! high cloud (less subgrid variability)
            rh0adj = rhmax-(rhmax-rh0(i,j))/(d_one+0.15D0*(tc0-t3(i,k,j)))
          end if
          if ( rh3(i,k,j) >= rhmax ) then        ! full cloud cover
            fcc(i,k,j) = d_one
          else if ( rh3(i,k,j) <= rh0adj ) then  ! no cloud cover
            fcc(i,k,j) = d_zero
          else                                   ! partial cloud cover
            fcc(i,k,j) = d_one-dsqrt(d_one-(rh3(i,k,j)-rh0adj) / &
                          (rhmax-rh0adj))
            fcc(i,k,j) = dmin1(dmax1(fcc(i,k,j),0.01D0),0.99D0)
          end if !  rh0 threshold
!---------------------------------------------------------------------
!         Correction:
!         Ivan Guettler, 14.10.2010.
!         Based on: Vavrus, S. and Waliser D., 2008, 
!         An Improved Parameterization for Simulating Arctic Cloud Amount
!         in the CCSM3 Climate Model, J. Climate 
!---------------------------------------------------------------------
          if ( p3(i,k,j) >= 75.0D0 ) then
            ! Clouds below 750hPa
            if ( qv3(i,k,j) <= 0.003D0 ) then
              fcc(i,k,j) = fcc(i,k,j) * &
                     dmax1(0.15D0,dmin1(d_one,qv3(i,k,j)/0.003D0))
            end if
          end if
!---------------------------------------------------------------------
!         End of the correction.
!---------------------------------------------------------------------
        end do
      end do

!--------------------------------------------------------------------
!     2.  Combine large-scale and convective fraction and liquid water
!         to be passed into radiation.
!--------------------------------------------------------------------
      do k = 1 , kz
        do i = 2 , iym2
          ! Cloud Water Volume
          ! kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
          exlwc = qc3(i,k,j)*rho3(i,k,j)*d_1000

          ! temperature dependance for convective cloud water content
          ! in g/m3 (Lemus et al., 1997)
          cldlwc(i,k)  = 0.127D+00 + 6.78D-03*(t3(i,k,j)-tzero) + &
                         1.29D-04* (t3(i,k,j)-tzero)**d_two  +    &
                         8.36D-07*(t3(i,k,j)-tzero)**d_three

          if ( cldlwc(i,k) > 0.3D+00 ) cldlwc(i,k) = 0.3D+00
          if ( (t3(i,k,j)-tzero) < -50D+00 ) cldlwc(i,k) = 0.001D+00
          ! Apply the parameterisation based on temperature to the
          ! convective fraction AND the large scale clouds :
          ! the large scale cloud water content is not really used by
          ! current radiation code, needs further evaluation.
          !TAO: but only apply this parameterization to large scale LWC 
          !if the user specifies it
          if (iconvlwp == 1) exlwc = cldlwc(i,k)
          cldlwc(i,k) = (cldfra(i,k)*cldlwc(i,k)+fcc(i,k,j)*exlwc) / &
                        dmax1(cldfra(i,k)+fcc(i,k,j),0.01D0)
          cldfra(i,k) = dmin1(dmax1(cldfra(i,k),fcc(i,k,j)),fcmax)
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
!    (see Pal et al 2000).  When partial clouds exist, the qvten  c
!    in/out of the clear and cloudy portions of the grid cell is  c
!    assumed to be at the same rate (i.e., if there is 80% cloud  c
!    cover, .2 of qvten goes to raising qv in the clear region    c
!    and .8 goes to condensation or evaporation of qc in the      c
!    cloudy portion).                                             c
!                                                                 c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
    subroutine condtq(jstart,jstop,psc,qvcs)
! 
    implicit none
!
    integer , intent(in) :: jstart , jstop
    real(8) , dimension(iy,jxp) , intent(in) :: psc
    real(8) , dimension(iy,kz) , intent(inout) :: qvcs
!
!   rhc    - Relative humidity at ktau+1
!   rh0adj - Adjusted relative humidity threshold at ktau+1
!   fccc   - Cloud fraction at ktau+1
!
    real(8) :: dqv , exces , fccc , pres , qvc_cld , qvs , &
               r1 , rh0adj , rhc , satvp
    integer :: i , j , k

    do j = jstart , jstop

!---------------------------------------------------------------------
!     1.  Compute t, qv, and qc at tau+1 without condensational term
!---------------------------------------------------------------------
      do k = 1 , kz
        do i = 2 , iym2
          tmp3(i,k) = (t2(i,k,j)+dt*tten(i,k,j))/psc(i,j)
          qvcs(i,k) = dmax1((qv2(i,k,j)+dt*qvten(i,k,j))/psc(i,j),lowq)
          qccs(i,k) = dmax1((qc2(i,k,j)+dt*qcten(i,k,j))/psc(i,j),lowq)
        end do
      end do
     
!---------------------------------------------------------------------
!     2.  Compute the cloud condensation/evaporation term.
!---------------------------------------------------------------------
      do k = 1 , kz
        do i = 2 , iym2
     
!         2a. Calculate the saturation mixing ratio and relative humidity
          pres = (a(k)*psc(i,j)+ptop)*d_1000
          if ( tmp3(i,k) > tzero ) then
            satvp = svp1*d_1000*dexp(svp2*(tmp3(i,k)-tzero)/(tmp3(i,k)-svp3))
          else
            satvp = svp4*d_1000*dexp(svp5-svp6/tmp3(i,k))
          end if
          qvs = dmax1(ep2*satvp/(pres-satvp),lowq)
          rhc = dmax1(qvcs(i,k)/qvs,lowq)

          r1 = d_one/(d_one+wlhv*wlhv*qvs/(rwat*cpd*tmp3(i,k)*tmp3(i,k)))
     
!         2b. Compute the relative humidity threshold at ktau+1
          if ( tmp3(i,k) > tc0 ) then
            rh0adj = rh0(i,j)
          else ! high cloud (less subgrid variability)
            rh0adj = rhmax - (rhmax-rh0(i,j))/(d_one+0.15D0*(tc0-tmp3(i,k)))
          end if
     
!         2c. Compute the water vapor in excess of saturation
          if ( rhc >= rhmax .or. rhc < rh0adj ) then ! Full or no cloud cover
            dqv = qvcs(i,k) - qvs*conf ! Water vapor in excess of sat
            tmp1(i,k) = r1*dqv
          else                                       ! Partial cloud cover
            fccc = d_one - dsqrt(d_one-(rhc-rh0adj)/(rhmax-rh0adj))
            fccc = dmin1(dmax1(fccc,0.01D0),d_one)
            qvc_cld = dmax1((qs3(i,k,j)+dt*qvten(i,k,j)/psc(i,j)),d_zero)
            dqv = qvc_cld - qvs*conf       ! qv diff between predicted qv_c
            tmp1(i,k) = r1*dqv*fccc        ! grid cell average
          end if
     
!         2d. Compute the new cloud water + old cloud water
          exces = qccs(i,k) + tmp1(i,k)
          if ( exces >= d_zero ) then ! Some cloud is left
            tmp2(i,k) = tmp1(i,k)/dt
          else                        ! The cloud evaporates
            tmp2(i,k) = -qccs(i,k)/dt
          end if
     
        end do
      end do
     
!---------------------------------------------------------------------
!     3.  Compute the tendencies.
!---------------------------------------------------------------------
      do k = 1 , kz
        do i = 2 , iym2
          qvten(i,k,j) = qvten(i,k,j) - psc(i,j)*tmp2(i,k)
          qcten(i,k,j) = qcten(i,k,j) + psc(i,j)*tmp2(i,k)
          tten(i,k,j) = tten(i,k,j) + psc(i,j)*tmp2(i,k)*wlhvocp
        end do
      end do
    end do
   
    end subroutine condtq
!
end module mod_precip
