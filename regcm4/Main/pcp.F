C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine pcp(j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     This subroutine computes the 'large scale' precipitation        c
c     based on the excess of cloud water above a threshold value.     c
c     The threshold (qcth) is temperature dependant and is based      c
c     on in cloud measurements of liquid cloud water (not ice).       c
c     Rain is only produced from the cloudy fraction of a grid cell   c
c     but the calculated precip value is a grid cell average.         c
c                                                                     c
c     This routine also computes raindrop evaporation and accretion.  c
c                                                                     c
c     See Pal et al. 2000 JGR-Atmos for more information.             c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer j
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/bats2.cb"
#     include "include/main.cb"
#     include "include/cvaria.cb"
#     include "include/pmoist.cb"
#     include "include/slice.cb"
#     include "include/trachem.cb"

      real(kind=8)  pptsum(ix)               ! Precip sum beginning from top
      integer i,k
      real(kind=8)  thog,dpovg,qcincld,pptnew,qcleft
      real(kind=8)  pptacc,rhcs,rdevap,uconv,aprdiv,i1000
      real(kind=8)  es,tk,tcel,q,qs,rh,afc,p,rho,qcw
      real(kind=8)  pptmax ! maximum precipation rate (total cloud water/dt)
      real(kind=8)  pptkm1 ! precipation accumulated from above
Cchem2
      integer kk
      real(kind=8)  uch
Cchem2_

c--------------------------------------------------------------------
c   1. Compute the precipitation formed in each layer.
c      The computations are performed from the top to the surface.
c      - Auto-conversion: similar to SIMEX (Giorgi and Shields 1999).
c      - Raindrop Accretion:  Beheng (1994); EQN 6
c      - Raindrop Evaporation:  Sundqvist (1988); EQN 3.4a
c--------------------------------------------------------------------


c   1a. Perform computations for the top layer (layer 1)
      thog = 1000./g
      i1000=1./1000.
Cchem2
      do k=1,kx
      do i=1,ix
        remrat(i,k) = 0.0
      end do
      end do
Cchem2_
      do i=2,ilxm

        afc = fcc(i,1,j)                                      ![frac][avg]

        if (afc.gt.0.01) then ! if there is a cloud
c     1aa. Compute temperature and humidities with the adjustments
c          due to convection.
c         q = qvb3d(i,1,j) + qcuten(i,1)*dt                  ![kg/kg][avg]
c         tk = tb3d(i,1,j) + tcuten(i,1)*dt                  ![k][avg]
          q = qvb3d(i,1,j)                                   ![kg/kg][avg]
          tk = tb3d(i,1,j)                                   ![k][avg]
          tcel = tk - 273.15                                 ![C][avg]
          p = pb3d(i,1,j)*1000.                              ![Pa][avg]
          rho=p/(r*tk)                                       ![kg/m3][avg]
          qcw = qcb3d(i,1,j)                                 ![kg/kg][avg]
c     1ab. Calculate the in cloud mixing ratio [kg/kg]
          qcincld = qcw/afc                                   ![kg/kg][cld]
c     1ac. Compute the maximum precipation rate
c          (i.e. total cloud water/dt) [kg/kg/s]
          pptmax = qcw/dt                                    ![kg/kg/s][avg]
c     1ad. Implement here the formula for qcth.
c          - Gultepe & Isaac, J. Clim, 1997, v10 p446 table 4, eq 5
c          - The factor of 1000 converts from g/kg to kg/kg
c          - The factor of cgul accounts for the fact that the Gultepe
c            and Isaac equation is for mean cloud water while qcth is the
c            theshhold for auto-conversion.
          qcth = cgul(i,j)*(10.**(-0.489+0.0134*tcel))*i1000 ![kg/kg][cld]
c     1ae. Compute the gridcell average autoconversion [kg/k g/s]
          pptnew = qck1(i,j)*(qcincld-qcth)*afc               ![kg/kg/s][avg]
          pptnew = dmin1(dmax1(pptnew,0.0d0),pptmax)         ![kg/kg/s][avg]
          if (pptnew.gt.0.0) then  ! New precipitation
c     1af. Compute the cloud removal rate (for chemistry) [1/s]
Cchem2
            remrat(i,1) = pptnew/qcw
Cchem2_
c     1ag. Compute the amount of cloud water removed by raindrop accretion
c          [kg/kg/s].  In the layer where the precipitation is formed,
c          only half of the precipitation is assumed to accrete.
c       1aga. Compute the amount of water remaining in the cloud [kg/kg]
            qcleft = qcw-pptnew*dt                           ![kg/kg][avg]
c       1agb. Add 1/2 of the new precipitation can accrete.
            pptkm1 = 0.5*pptnew/afc*rho*dt                    ![kg/m3][cld]
c       1agc. Accretion [kg/kg/s]=[m3/kg/s]*[kg/kg]*[kg/m3]
            pptacc = caccr*qcleft*pptkm1                     ![kg/kg/s][avg]
c       1agd. Update the precipitation accounting for the accretion [kg/kg/s]
            pptnew = dmin1(pptmax,pptacc+pptnew)             ![kg/kg/s][avg]
c     1ah. Accumulate precipitation and convert to kg/m2/s
            dpovg = dsigma(1)*psb(i,j)*thog                  ![kg/m2]
            pptsum(i) = pptnew*dpovg                         ![kg/m2/s][avg]
c     1ai. Compute the cloud water tendency [kg/kg/s*cb]
            qcten(i,1,j) = qcten(i,1,j) - pptnew*psb(i,j)    ![kg/kg/s*cb][avg]
          else  ! Cloud but no new precipitation
            pptsum(i) = 0.0                                  ![kg/m2/s][avg]
          end if
        else  ! No cloud
          pptsum(i) = 0.0                                    ![kg/m2/s][avg]
        end if

      end do

C ****  LAYER TWO TO KL ****
c   1b. Perform computations for the 2nd layer to the surface
      do k=2,kx
      do i=2,ilxm

c     1ba. Compute temperature and humidities with the adjustments
c          due to convection.
c       q = qvb3d(i,k,j) ! + qcuten(i,k)*dt                  ![kg/kg][avg]
c       tk = tb3d(i,k,j) ! + tcuten(i,k)*dt                  ![k][avg]
        q = qvb3d(i,k,j)                                     ![kg/kg][avg]
        tk = tb3d(i,k,j)                                     ![k][avg]
        tcel = tk - 273.15                                   ![C][avg]
        p = pb3d(i,k,j)*1000.                                ![Pa][avg]
        rho=p/(r*tk)                                         ![kg/m3][avg]
        qcw = qcb3d(i,k,j)                                   ![kg/kg][avg]
        afc = fcc(i,k,j)                                      ![frac][avg]
        if (tcel.gt.0.0) then
          es=svp1*1000.*dexp(svp2*tcel/(tk-svp3))            ![Pa][avg]
        else
          es=svp1*1000.*dexp(22.514-6.15e3/tk)               ![Pa][avg]
        end if
        qs=ep2*es/(p-es)                                     ![kg/kg][avg]
        rh = dmin1(dmax1(q/qs,0.0d0),rhmax)                  ![frac][avg]

c     1bb. Convert accumlated precipitation to kg/kg/s.
c          Used for raindrop evaporation and accretion.
        dpovg = dsigma(k)*psb(i,j)*thog                      ![kg/m2][avg]
        pptkm1 = pptsum(i)/dpovg                             ![kg/kg/s][avg]

c     1bc. Compute the raindrop evaporation in the clear portion of
c          the gridcell.
c          - It is assumed that raindrops do not evaporate in clouds
c            and the rainfall from above is evenly distributed in
c            gridcell (i.e. the gridcell average precipitation is used).
        if (pptsum(i).gt.0.0 .and. afc.lt.0.99) then
c       2bca. Compute the clear sky relative humidity
          rhcs = (rh-afc*rhmax)/(1.0-afc)                      ![frac][clr]
          rhcs = dmax1(dmin1(rhcs,rhmax),0.0d0)              ![frac][clr]
c       2bcb. Raindrop evaporation [kg/kg/s]
          rdevap = cevap*(rhmax-rhcs)*sqrt(pptsum(i))*(1.-afc)![kg/kg/s][avg]
          rdevap = dmin1((qs-q)/dt,rdevap)                   ![kg/kg/s][avg]
          rdevap = dmin1(dmax1(rdevap,0.0d0),pptkm1)         ![kg/kg/s][avg]
c       2bcc. Update the precipitation accounting for the raindrop
c          evaporation [kg/m2/s]
          pptsum(i) = pptsum(i)-rdevap*dpovg                 ![kg/m2/s][avg]
c       2bcf. Compute the water vapor tendency [kg/kg/s*cb]
          qvten(i,k,j) = qvten(i,k,j) + rdevap*psb(i,j)      ![kg/kg/s*cb][avg]
c       2bcf. Compute the temperature tendency [K/s*cb]
          tten(i,k,j) = tten(i,k,j) - xlvocp*rdevap*psb(i,j) ![k/s*cb][avg]
        else  ! no precipitation from above
          rdevap = 0.0                                       ![kg/kg/s][avg]
        end if

c     1bd. Compute the autoconversion and accretion [kg/kg/s]
        if (afc.gt.0.01) then ! if there is a cloud
c       1bda. Calculate the in cloud mixing ratio [kg/kg]
          qcincld = qcw/afc                                   ![kg/kg][cld]
c       1bdb. Compute the maximum precipation rate
c             (i.e. total cloud water/dt) [kg/kg/s]
          pptmax = qcw/dt                                    ![kg/kg/s][cld]
c       1bdc. Implement the Gultepe & Isaac formula for qcth.
          qcth = cgul(i,j)*(10.**(-0.489+0.0134*tcel))*i1000 ![kg/kg][cld]
c       1bdd. Compute the gridcell average autoconversion [kg/kg/s]
          pptnew = qck1(i,j)*(qcincld-qcth)*afc               ![kg/kg/s][avg]
          pptnew = dmin1(dmax1(pptnew,0.0d0),pptmax)         ![kg/kg/s][avg]
c     1be. Compute the cloud removal rate (for chemistry) [1/s]
Cchem2
          if (pptnew.gt.0.0) remrat(i,k) = pptnew/qcw
Cchem2_

c     1bf. Compute the amount of cloud water removed by raindrop accretion
c          [kg/kg/s].  In the layer where the precipitation is formed,
c          only half of the precipitation can accrete.
          if (pptkm1.gt.0.0 .or. pptnew.gt.0.0) then
c       1bfa. Compute the amount of water remaining in the cloud [kg/kg]
            qcleft = dmax1(qcw-pptnew*dt,0.d0)               ![kg/kg][avg]
c       1bfb. Add 1/2 of the new precipitation to the accumulated
c             precipitation [kg/m3]
            pptkm1 = (pptkm1+0.5*pptnew/afc)*rho*dt           ![kg/m3][cld]
c       1bfc. accretion [kg/kg/s]
            pptacc = caccr*qcleft*pptkm1                     ![kg/kg/s][avg]
c       1bfd. Update the precipitation accounting for the accretion [kg/kg/s]
            pptnew = dmin1(pptmax,pptacc+pptnew)             ![kg/kg/s][avg]
          end if
c     1bg. Accumulate precipitation and convert to kg/m2/s
          pptsum(i) = pptsum(i) + pptnew*dpovg               ![kg/m2/s][avg]
c     1bh. Compute the cloud water tendency [kg/kg/s*cb]
          qcten(i,k,j) = qcten(i,k,j) - pptnew*psb(i,j)      ![kg/kg/s*cb][avg]
        else
          pptnew = 0.0                                       ![kg/kg/s][avg]
        end if

      end do
      end do

Cchem2
c--------------------------------------------------------------------
c   2. Perform aerosol removal computations
c      - swith do i,k loop, add rembc (the below cloud scavenging
c        rate, s^-1)
c      - Levin & Schwatz
c--------------------------------------------------------------------
      if (ichem.eq.1) then
        uch = 1000./g*3600.
        do 33 i=2,ilxm
          rembc(i,1)=0.
          do 34 k=2,kx
            rembc(i,k)=0.
            if(remrat(i,k).gt.0.) then
              do kk=1,k-1
                rembc(i,k) = rembc(i,k)         ! mm/hr
     &                     + remrat(i,kk)*qcb3d(i,kk,j)*psb(i,j)
     &                     * dsigma(kk)*uch
              end do
              rembc(i,k)=6.5*1.e-5*rembc(i,k)**.68     ! s^-1
            end if
34        continue
33      continue
      end if
Cchem2_


c--------------------------------------------------------------------
c   3. Convert the accumlated precipitation to appropriate units for
c      the surface physics and the output
c--------------------------------------------------------------------
      uconv = 60.*dtmin
      aprdiv = 1./dble(nbatst)
      if (jyear.eq.jyear0.and.ktau.eq.0) aprdiv = 1.
      do 60 i=2,ilxm
        rainnc(i,j) = rainnc(i,j) + pptsum(i)*uconv
        pptnc(i,j) = pptnc(i,j) + pptsum(i)*aprdiv
 60   continue

      return
      end subroutine pcp
