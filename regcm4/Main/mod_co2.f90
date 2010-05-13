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
 
      module mod_co2

      implicit none

      contains

      subroutine co2

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     model for plant carbon uptake and dead carbon decomposition
!
!            dtbat = surface output interval
!            cari = rate of carbon uptake by plants
!       pbp1d(n,i) = accumulated primary biomass
!      resp1d(n,i) = carbon in soil and in veg (kg c / m**2 / sec)
!              ra = leaf aerodynamic resistance factor
!             rap = leaf boundary layer resistance to co2
!           resps = soil respiration
!             rmp = physical mesophyllic resistance
!              rs = stomatal resistance
!             rsp = stomatal resistance to co2
!              rt = total mechanical resistance to co2
!
!     0.33 pptv co2 = 0.30 pptm dry matter in plants
!
!     decomposable soil carbon assumed incremented by npp-soil resp.
!
!     resistances for co2 are larger than those for h2o due to
!                    difference in molecular weight
!
      use mod_dynparam
      use mod_param1 , only : dtbat
      use mod_bats , only : fdry , pbp1d , resp1d , ldoc1d , sigf ,     &
                   & rlai , gwmx0, solis , tlef1d , ssw1d , tg1d ,      &
                   & xlai , xlsai
      use mod_ictp01
      implicit none
!
! Local variables
!
      real(8) , dimension(nnsg,iym1) :: cari
      integer :: n , i
      real(8) :: rap , resps , rmp , rsp , rt
! 
      rmp = 800.0
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              rsp = rs(n,i)*1.7
              rap = ra(n,i)*1.5
              rt = rsp + rap + rmp
              cari(n,i) = sigf(n,i)*xlsai(n,i)*fdry(n,i)                &
                         & *carbon(solis(i)*rlai(n,i),tlef1d(n,i),rt,   &
                         & tg1d(n,i),xlai(n,i),xlsai(n,i))
              pbp1d(n,i) = pbp1d(n,i) + cari(n,i)*dtbat
            end if
          end if
        end do
      end do
 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              if ( pbp1d(n,i).lt.0 ) pbp1d(n,i) = 0.
              resps = 0.7E-7*resp1d(n,i)*dexp(0.1*(tg1d(n,i)-300.))     &
                    & *dmin1(1.D0,ssw1d(n,i)/(0.6*gwmx0(n,i)))
              resp1d(n,i) = resp1d(n,i) + (cari(n,i)-resps)*dtbat
              if ( resp1d(n,i).lt.0.0 ) resp1d(n,i) = 0.
            end if
          end if
        end do
      end do
      end subroutine co2
!
!
!
      function carbon(vf,t,rm,tg,xlai,xlsai)
      implicit none
!
! Dummy arguments
!
      real(8) :: rm , t , tg , vf , xlai , xlsai
      real(8) :: carbon
      intent (in) rm , tg , vf , xlai , xlsai
!
! Local variables
!
      real(8) :: a , ab , ac , al , alphtl , b , bc , betatl , cco2 ,   &
               & cco2i , ccold , gt , p , pm , pml , rt , sl , tmx , w ,&
               & wd , wp , xk , xkb , xl
      real(8) :: e , g , r
      integer :: it
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ***  a model for whole leaf photosynthesis based on that derived
!     *** by tenhunen as summarized in the text of d gates
!
!     function g represents temperature dependence of photosynthesis
!======================================================================
!     still to improve:
!     parameter for different vegetation types, now wheat only
!     correct initial values for resp must be set in initb
!     parameterisation of resps for cultivated land ok ?
!     betatl = rai/lai  (must be set reasonably)
!====================================================================
!     this fn is not in si - vf (radn) and rm (res) handed in are in si
!     result of fn (carbon) is handed back in si (kg/m**2/s)
!======================================================================
      g(t,tmx,sl) = dexp(sl*(1./tmx-1./t))                              &
                  & /(1.+(dexp(sl*(1./tmx-1./t)*6)))*5.E-3*t
!     ***  temperature dependence of dark respiration
      r(t) = dexp(30.-9.E3/t)
!     ****  light dependence of photosynthesis
      e(xl,a,pml) = a*xl/(1.+(a*xl/pml)**2)**0.5
 
      p = 0
      rt = r(t)
!     ****    alphtl=sai/lai     betatl=rai/lai  (rai=root area index)
      alphtl = (xlsai-xlai)/xlai
      betatl = 0.5
      if ( vf.lt.2. ) then
 
!       ***   nighttime maintenance respiration only
        carbon = -0.36E-8*((1.+alphtl)*0.877*rt+betatl*(r(tg)-0.123*rt))
      else
!       **** convert lambda less than 0.7 micron solar into photon units
!       *****  light intensity (e/m2)
        xl = 4.6E-3*vf
!       *** co2 external concentration(mm/m3)
!       ****    335ppm/v
        cco2 = 13.5
!       ***initial guess for co2 concentration inside chloroplast
        cco2i = cco2
!       ****  co2 half max in absence of oxygen  (mm/m3)
        xk = 0.5
!       ****  oxygen inhibition factor
        b = 3.56
        xkb = xk*b
!       ****  maximum temperature optimum light saturated photosynthesis
        pml = 0.050
!       ****    quantum yield(mm/e)
        al = 0.05
        gt = g(t,320.D0,4.D3)
!       ****  maximum photosynthesis
        pm = e(xl,al,pml*gt)
!       ******   iterate
 
        do it = 1 , 30
!         ****  photorespiration
          wp = pm/(1.+0.4*(1.+cco2i/xk))
!         ****    total respiration
!         ****dark respiration within daytime leaves
          wd = 3.E-4*rt + 0.14*p
          w = wp + wd
!         *****    carbon uptake factors
          ac = cco2 + xkb + rm*(pm-w)
          bc = 4.*rm*(cco2*(pm-w)-xkb*w)
          ab = dsqrt(ac**2-bc)
          p = 0.5*(ac-ab)/rm
          ccold = cco2i
          cco2i = cco2 - p*rm
          if ( dabs(cco2i-ccold).le.0.05 .and. it.gt.9 ) exit
        end do
!       **** respiration outside leaves
        carbon = 1.2E-5*((1.-0.14*(alphtl+betatl))                      &
               & *p-3.0E-4*(alphtl*rt+betatl*r(tg)))
        return
      end if
 
      end function carbon
!
      end module mod_co2
