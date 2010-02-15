!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine pcp(j)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     This subroutine computes the 'large scale' precipitation        c
!     based on the excess of cloud water above a threshold value.     c
!     The threshold (qcth) is temperature dependant and is based      c
!     on in cloud measurements of liquid cloud water (not ice).       c
!     Rain is only produced from the cloudy fraction of a grid cell   c
!     but the calculated precip value is a grid cell average.         c
!                                                                     c
!     This routine also computes raindrop evaporation and accretion.  c
!                                                                     c
!     See Pal et al. 2000 JGR-Atmos for more information.             c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_param3 , only : dsigma
      use mod_bats
      use mod_main
      use mod_cvaria
      use mod_pmoist
      use mod_slice
      use mod_trachem
      use mod_constants , only : rgti , rgas , ep2 , wlhvocp , svp1 ,   &
                               & svp2 , svp3 , tmelt
      implicit none
!
! Dummy arguments
!
      integer :: j
      intent (in) j
!
! Local variables
!
      real(8) :: aprdiv , dpovg , es , afc , i1000 , p , pptacc ,       &
               & pptkm1 , pptmax , pptnew , q , qcincld , qcleft , qcw ,&
               & qs , rdevap , rh , rhcs , rho , tcel , thog , tk ,     &
               & uch , uconv
      integer :: i , k , kk
      real(8) , dimension(ix) :: pptsum
!
 
! Precip sum beginning from top
! maximum precipation rate (total cloud water/dt)
 
!--------------------------------------------------------------------
!     1. Compute the precipitation formed in each layer.
!     The computations are performed from the top to the surface.
!     - Auto-conversion: similar to SIMEX (Giorgi and Shields 1999).
!     - Raindrop Accretion:  Beheng (1994); EQN 6
!     - Raindrop Evaporation:  Sundqvist (1988); EQN 3.4a
!--------------------------------------------------------------------
 
 
!     1a. Perform computations for the top layer (layer 1)
      thog = 1000.*rgti       ! precipation accumulated from above
      i1000 = 1./1000.
!chem2
      do k = 1 , kx
        do i = 1 , ix
          remrat(i,k) = 0.0
        end do
      end do
!chem2_
      do i = 2 , ixm2
 
        afc = fcc(i,1,j)                                      ![frac][avg]
 
        if ( afc.gt.0.01 ) then
                             ! if there is a cloud
!         1aa. Compute temperature and humidities with the adjustments
!         due to convection.
!         q = qvb3d(i,1,j) + qcuten(i,1)*dt                 
!         ![kg/kg][avg] tk = tb3d(i,1,j) + tcuten(i,1)*dt              
!         ![k][avg]
          q = qvb3d(i,1,j)                                   ![kg/kg][avg]
          tk = tb3d(i,1,j)                                   ![k][avg]
          tcel = tk - tmelt                                  ![C][avg]
          p = pb3d(i,1,j)*1000.                              ![Pa][avg]
          rho = p/(rgas*tk)                                  ![kg/m3][avg]
          qcw = qcb3d(i,1,j)                                 ![kg/kg][avg]
!         1ab. Calculate the in cloud mixing ratio [kg/kg]
          qcincld = qcw/afc                                  ![kg/kg][cld]
!         1ac. Compute the maximum precipation rate
!         (i.e. total cloud water/dt) [kg/kg/s]
          pptmax = qcw/dt                                    ![kg/kg/s][avg]
!         1ad. Implement here the formula for qcth.
!         - Gultepe & Isaac, J. Clim, 1997, v10 p446 table 4, eq 5
!         - The factor of 1000 converts from g/kg to kg/kg
!         - The factor of cgul accounts for the fact that the Gultepe
!         and Isaac equation is for mean cloud water while qcth is the
!         theshhold for auto-conversion.
          qcth = cgul(i,j)*(10.**(-0.489+0.0134*tcel))*i1000 ![kg/kg][cld]
!         1ae. Compute the gridcell average autoconversion [kg/k g/s]
          pptnew = qck1(i,j)*(qcincld-qcth)*afc              ![kg/kg/s][avg]
          pptnew = dmin1(dmax1(pptnew,0.0D0),pptmax)         ![kg/kg/s][avg]
          if ( pptnew.gt.0.0 ) then
                                   ! New precipitation
!           1af. Compute the cloud removal rate (for chemistry) [1/s]
!chem2
            remrat(i,1) = pptnew/qcw
!chem2_
!           1ag. Compute the amount of cloud water removed by raindrop
!           accretion [kg/kg/s].  In the layer where the precipitation
!           is formed, only half of the precipitation is assumed to
!           accrete. 1aga. Compute the amount of water remaining in the
!           cloud [kg/kg]
            qcleft = qcw - pptnew*dt                         ![kg/kg][avg]
!           1agb. Add 1/2 of the new precipitation can accrete.
            pptkm1 = 0.5*pptnew/afc*rho*dt                   ![kg/m3][cld]
!           1agc. Accretion [kg/kg/s]=[m3/kg/s]*[kg/kg]*[kg/m3]
            pptacc = caccr*qcleft*pptkm1                     ![kg/kg/s][avg]
!           1agd. Update the precipitation accounting for the accretion
!           [kg/kg/s]
            pptnew = dmin1(pptmax,pptacc+pptnew)             ![kg/kg/s][avg]
!           1ah. Accumulate precipitation and convert to kg/m2/s
            dpovg = dsigma(1)*psb(i,j)*thog                  ![kg/m2]
            pptsum(i) = pptnew*dpovg                         ![kg/m2/s][avg]
!           1ai. Compute the cloud water tendency [kg/kg/s*cb]
            qcten(i,1,j) = qcten(i,1,j) - pptnew*psb(i,j)    ![kg/kg/s*cb][avg]
          else  ! Cloud but no new precipitation
            pptsum(i) = 0.0                                  ![kg/m2/s][avg]
          end if
        else  ! No cloud
          pptsum(i) = 0.0                                    ![kg/m2/s][avg]
        end if
 
      end do
 
!     ****  LAYER TWO TO KL ****
!     1b. Perform computations for the 2nd layer to the surface
      do k = 2 , kx
        do i = 2 , ixm2
 
!         1ba. Compute temperature and humidities with the adjustments
!         due to convection.
!         q = qvb3d(i,k,j) ! + qcuten(i,k)*dt                 
!         ![kg/kg][avg] tk = tb3d(i,k,j) ! + tcuten(i,k)*dt            
!         ![k][avg]
          q = qvb3d(i,k,j)                                   ![kg/kg][avg]
          tk = tb3d(i,k,j)                                   ![k][avg]
          tcel = tk - tmelt                                  ![C][avg]
          p = pb3d(i,k,j)*1000.                              ![Pa][avg]
          rho = p/(rgas*tk)                                  ![kg/m3][avg]
          qcw = qcb3d(i,k,j)                                 ![kg/kg][avg]
          afc = fcc(i,k,j)                                   ![frac][avg]
          if ( tcel.gt.0.0 ) then
            es = svp1*1000.*dexp(svp2*tcel/(tk-svp3))        ![Pa][avg]
          else
            es = svp1*1000.*dexp(22.514-6.15E3/tk)           ![Pa][avg]
          end if
          qs = ep2*es/(p-es)                                 ![kg/kg][avg]
          rh = dmin1(dmax1(q/qs,0.0D0),rhmax)                ![frac][avg]
 
!         1bb. Convert accumlated precipitation to kg/kg/s.
!         Used for raindrop evaporation and accretion.
          dpovg = dsigma(k)*psb(i,j)*thog                    ![kg/m2][avg]
          pptkm1 = pptsum(i)/dpovg                           ![kg/kg/s][avg]
 
!         1bc. Compute the raindrop evaporation in the clear portion of
!         the gridcell.
!         - It is assumed that raindrops do not evaporate in clouds
!         and the rainfall from above is evenly distributed in
!         gridcell (i.e. the gridcell average precipitation is used).
          if ( pptsum(i).gt.0.0 .and. afc.lt.0.99 ) then
!           2bca. Compute the clear sky relative humidity
            rhcs = (rh-afc*rhmax)/(1.0-afc)                  ![frac][clr]
            rhcs = dmax1(dmin1(rhcs,rhmax),0.0D0)            ![frac][clr]
!           2bcb. Raindrop evaporation [kg/kg/s]
            rdevap = cevap*(rhmax-rhcs)*sqrt(pptsum(i))*(1.-afc)
                                                             ![kg/kg/s][avg]
            rdevap = dmin1((qs-q)/dt,rdevap)                 ![kg/kg/s][avg]
            rdevap = dmin1(dmax1(rdevap,0.0D0),pptkm1)       ![kg/kg/s][avg]
!           2bcc. Update the precipitation accounting for the raindrop
!           evaporation [kg/m2/s]
            pptsum(i) = pptsum(i) - rdevap*dpovg             ![kg/m2/s][avg]
!           2bcf. Compute the water vapor tendency [kg/kg/s*cb]
            qvten(i,k,j) = qvten(i,k,j) + rdevap*psb(i,j)    ![kg/kg/s*cb][avg]
!           2bcf. Compute the temperature tendency [K/s*cb]
            tten(i,k,j) = tten(i,k,j) - wlhvocp*rdevap*psb(i,j)
                                                             ![k/s*cb][avg]
          else
              ! no precipitation from above
            rdevap = 0.0                                     ![kg/kg/s][avg]
          end if
 
!         1bd. Compute the autoconversion and accretion [kg/kg/s]
          if ( afc.gt.0.01 ) then
                             ! if there is a cloud
!           1bda. Calculate the in cloud mixing ratio [kg/kg]
            qcincld = qcw/afc                                ![kg/kg][cld]
!           1bdb. Compute the maximum precipation rate
!           (i.e. total cloud water/dt) [kg/kg/s]
            pptmax = qcw/dt                                  ![kg/kg/s][cld]
!           1bdc. Implement the Gultepe & Isaac formula for qcth.
            qcth = cgul(i,j)*(10.**(-0.489+0.0134*tcel))*i1000
                                                             ![kg/kg][cld]
!           1bdd. Compute the gridcell average autoconversion [kg/kg/s]
            pptnew = qck1(i,j)*(qcincld-qcth)*afc            ![kg/kg/s][avg]
            pptnew = dmin1(dmax1(pptnew,0.0D0),pptmax)       ![kg/kg/s][avg]
!           1be. Compute the cloud removal rate (for chemistry) [1/s]
!chem2
            if ( pptnew.gt.0.0 ) remrat(i,k) = pptnew/qcw
!chem2_
 
!           1bf. Compute the amount of cloud water removed by raindrop
!           accretion [kg/kg/s].  In the layer where the precipitation
!           is formed, only half of the precipitation can accrete.
            if ( pptkm1.gt.0.0 .or. pptnew.gt.0.0 ) then
!             1bfa. Compute the amount of water remaining in the cloud
!             [kg/kg]
              qcleft = dmax1(qcw-pptnew*dt,0.D0)             ![kg/kg][avg]
!             1bfb. Add 1/2 of the new precipitation to the accumulated
!             precipitation [kg/m3]
              pptkm1 = (pptkm1+0.5*pptnew/afc)*rho*dt        ![kg/m3][cld]
!             1bfc. accretion [kg/kg/s]
              pptacc = caccr*qcleft*pptkm1                   ![kg/kg/s][avg]
!             1bfd. Update the precipitation accounting for the
!             accretion [kg/kg/s]
              pptnew = dmin1(pptmax,pptacc+pptnew)           ![kg/kg/s][avg]
            end if
!           1bg. Accumulate precipitation and convert to kg/m2/s
            pptsum(i) = pptsum(i) + pptnew*dpovg             ![kg/m2/s][avg]
!           1bh. Compute the cloud water tendency [kg/kg/s*cb]
            qcten(i,k,j) = qcten(i,k,j) - pptnew*psb(i,j)    ![kg/kg/s*cb][avg]
          else
            pptnew = 0.0                                     ![kg/kg/s][avg]
          end if
 
        end do
      end do
 
!chem2
!--------------------------------------------------------------------
!     2. Perform aerosol removal computations
!     - swith do i,k loop, add rembc (the below cloud scavenging
!     rate, s^-1)
!     - Levin & Schwatz
!--------------------------------------------------------------------
      if ( ichem.eq.1 ) then
        uch = 1000.*rgti*3600.
        do i = 2 , ixm2
          rembc(i,1) = 0.
          do k = 2 , kx
            rembc(i,k) = 0.
            if ( remrat(i,k).gt.0. ) then
              do kk = 1 , k - 1
                rembc(i,k) = rembc(i,k) + remrat(i,kk)*qcb3d(i,kk,j)    &
                           & *psb(i,j)*dsigma(kk)*uch
                                                ! mm/hr
              end do
              rembc(i,k) = 6.5*1.E-5*rembc(i,k)**.68   ! s^-1
            end if
          end do
        end do
      end if
!chem2_
 
 
!--------------------------------------------------------------------
!     3. Convert the accumlated precipitation to appropriate units for
!     the surface physics and the output
!--------------------------------------------------------------------
      uconv = 60.*dtmin
      aprdiv = 1./dble(nbatst)
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) aprdiv = 1.
      do i = 2 , ixm2
        rainnc(i,j) = rainnc(i,j) + pptsum(i)*uconv
        pptnc(i,j) = pptnc(i,j) + pptsum(i)*aprdiv
      end do
 
      end subroutine pcp
