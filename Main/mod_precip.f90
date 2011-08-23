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
! Large Scale Precipitation -- Pal et al. 2000 JGR-Atmos
!
  use mod_runparams
  use mod_atm_interface , only : atmstate , slice , surfpstate , surfstate
!
  private
!
! Precip sum beginning from top
  real(8) , pointer , dimension(:) :: pptsum
  real(8) , pointer , dimension(:,:) :: psf , rainnc , lsmrnc
  real(8) , pointer , dimension(:,:,:) :: t3 , p3 , qv3 , qc3
  real(8) , pointer , dimension(:,:,:) :: tten , qvten , qcten
  real(8) , pointer , dimension(:,:) :: chrmrat , chrmbc

  real(8) :: qcth

  integer :: istart , istopx
!
  real(8) , parameter :: uch = d_1000*regrav*secph
!
  real(8) , public , pointer , dimension(:,:,:) :: fcc
  real(8) , public , pointer , dimension(:,:) :: qck1 , cgul
  real(8) , public :: caccr , cevap , rhmax
!
  public :: allocate_mod_precip , init_precip , pcp
!
  contains
!
    subroutine allocate_mod_precip
      implicit none
      call getmem3d(fcc,1,iy,1,kz,1,jxp,'pcp:fcc')
      call getmem2d(qck1,1,iy,1,jxp,'pcp:qck1')
      call getmem2d(cgul,1,iy,1,jxp,'pcp:cgul')
    end subroutine allocate_mod_precip
!
    subroutine init_precip(atmslice,tendency,sps,surface,lsmrainnc, &
                           chremrat,chrembc)
      implicit none
      type(slice) , intent(in) :: atmslice
      type(atmstate) , intent(in) :: tendency
      type(surfpstate) , intent(in) :: sps
      type(surfstate) , intent(in) :: surface
      real(8) , pointer , dimension(:,:) :: lsmrainnc
      real(8) , pointer , dimension(:,:) :: chremrat , chrembc

      t3      => atmslice%tb3d
      p3      => atmslice%pb3d
      qv3     => atmslice%qvb3d
      qc3     => atmslice%qcb3d
      tten    => tendency%t
      qcten   => tendency%qc
      qvten   => tendency%qv
      psf     => sps%ps
      rainnc  => surface%rainnc
      lsmrnc  => lsmrainnc
      chrmrat => chremrat
      chrmbc  => chrembc

      istart = lbound(rainnc,1) + 1
      istopx = ubound(rainnc,1) - 2

      call getmem1d(pptsum,istart,istopx,'pcp:pptsum')
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
    chrmrat(istart:istopx,1:kz) = d_zero

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
            chrmrat(i,1) = pptnew/qcw
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
            if ( pptnew > d_zero ) chrmrat(i,k) = pptnew/qcw
     
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
     
!--------------------------------------------------------------------
!     2. Perform aerosol removal computations
!       - swith do i,k loop, add chrmbc (the below cloud scavenging rate, s^-1)
!       - Levin & Schwatz
!--------------------------------------------------------------------
      if ( ichem == 1 ) then
        do i = istart , istopx
          chrmbc(i,1) = d_zero
          do k = 2 , kz
            chrmbc(i,k) = d_zero
            if ( chrmrat(i,k) > d_zero ) then
              do kk = 1 , k - 1
                chrmbc(i,k) = chrmbc(i,k) + chrmrat(i,kk)*qc3(i,kk,j) *  &
                             psf(i,j)*dsigma(kk)*uch          ![mm/hr]
              end do
              chrmbc(i,k) = 6.5D0*1.0D-5*chrmbc(i,k)**0.68D0   ![s^-1]
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
end module mod_precip
