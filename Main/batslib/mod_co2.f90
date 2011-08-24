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
!
  use mod_dynparam
  use mod_bats
  use mod_leaftemp
!
  private
!
  public :: co2
!
  contains
!
!=======================================================================
!l  based on: bats version 1e          copyright 18 august 1989
!=======================================================================
!
! CO2
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
!=======================================================================
!
  subroutine co2
!
    implicit none
!
    real(8) , dimension(nnsg,iym1) :: cari
    integer :: n , i
    real(8) :: rap , resps , rsp , rt , rcar
    real(8) , parameter :: rmp = 800.0D0
!
    do i = 2 , iym1
      do n = 1 , nnsg
        if ( ldoc1d(n,i) /= 0 ) then
          if ( sigf(n,i) > 0.001D0 ) then
            rsp = lftrs(n,i)*1.7D0
            rap = lftra(n,i)*1.5D0
            rt = rsp + rap + rmp
            rcar = carbon(solis(i)*rlai(n,i),tlef1d(n,i),rt, &
                          tg1d(n,i),xlai(n,i),xlsai(n,i))
            cari(n,i) = sigf(n,i)*xlsai(n,i)*fdry(n,i)*rcar
            pbp1d(n,i) = pbp1d(n,i) + cari(n,i)*dtbat
          end if
        end if
      end do
    end do
   
    do i = 2 , iym1
      do n = 1 , nnsg
        if ( ldoc1d(n,i) /= 0 ) then
          if ( sigf(n,i) > 0.001D0 ) then
            if ( pbp1d(n,i) < 0 ) pbp1d(n,i) = d_zero
            resps = 0.7D-7*resp1d(n,i)*dexp(0.1D0*(tg1d(n,i)-300.0D0)) * &
                    dmin1(d_one,ssw1d(n,i)/(0.6D0*gwmx0(n,i)))
            resp1d(n,i) = resp1d(n,i) + (cari(n,i)-resps)*dtbat
            if ( resp1d(n,i) < d_zero ) resp1d(n,i) = d_zero
          end if
        end if
      end do
    end do
  end subroutine co2
!
!=======================================================================
! CARBON
!
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
!
  function carbon(vf,t,rm,tg,xlai,xlsai)
    implicit none
!
    real(8) :: rm , t , tg , vf , xlai , xlsai
    real(8) :: carbon
    intent (in) rm , tg , vf , xlai , xlsai
!
    real(8) :: a , ab , ac , al , alphtl , b , bc , betatl , cco2 ,   &
               cco2i , ccold , gt , p , pm , pml , rt , sl , tmx , w ,&
               wd , wp , xk , xkb , xl
    real(8) :: e , g , r
    integer :: it
!
!====================================================================
!     this fn is not in si - vf (radn) and rm (res) handed in are in si
!     result of fn (carbon) is handed back in si (kg/m**2/s)
!======================================================================
!
    g(t,tmx,sl) = dexp(sl*(d_one/tmx-d_one/t)) / &
                 (d_one+(dexp(sl*(d_one/tmx-d_one/t)*6.0D0)))*5.0D-3*t
!     ***  temperature dependence of dark respiration
    r(t) = dexp(30.0D0-9.0D3/t)
!     ****  light dependence of photosynthesis
    e(xl,a,pml) = a*xl/(d_one+(a*xl/pml)**d_two)**d_half
   
    p = d_zero
    rt = r(t)
!     ****    alphtl=sai/lai     betatl=rai/lai  (rai=root area index)
    alphtl = (xlsai-xlai)/xlai
    betatl = d_half
    if ( vf < d_two ) then
   
!       ***   nighttime maintenance respiration only
      carbon = -0.36D-8*((d_one+alphtl)*0.877D0*rt+ &
                 betatl*(r(tg)-0.123D0*rt))
    else
!       **** convert lambda less than 0.7 micron solar into photon units
!       *****  light intensity (e/m2)
      xl = 4.6D-3*vf
!       *** co2 external concentration(mm/m3)
!       ****    335ppm/v
      cco2 = 13.5D0
!       ***initial guess for co2 concentration inside chloroplast
      cco2i = cco2
!       ****  co2 half max in absence of oxygen  (mm/m3)
      xk = d_half
!       ****  oxygen inhibition factor
      b = 3.56D0
      xkb = xk*b
!       ****  maximum temperature optimum light saturated photosynthesis
      pml = 0.050D0
!       ****    quantum yield(mm/e)
      al = 0.05D0
      gt = g(t,320.0D0,4.0D3)
!       ****  maximum photosynthesis
      pm = e(xl,al,pml*gt)
!       ******   iterate
   
      do it = 1 , 30
!         ****  photorespiration
        wp = pm/(d_one+0.4D0*(d_one+cco2i/xk))
!         ****    total respiration
!         ****dark respiration within daytime leaves
        wd = 3.0D-4*rt + 0.14D0*p
        w = wp + wd
!         *****    carbon uptake factors
        ac = cco2 + xkb + rm*(pm-w)
        bc = d_four*rm*(cco2*(pm-w)-xkb*w)
        ab = dsqrt(ac**d_two-bc)
        p = d_half*(ac-ab)/rm
        ccold = cco2i
        cco2i = cco2 - p*rm
        if ( dabs(cco2i-ccold) <= 0.05D0 .and. it > 9 ) exit
      end do
!       **** respiration outside leaves
      carbon = 1.2D-5*((d_one-0.14D0*(alphtl+betatl)) * &
               p-3.0D-4*(alphtl*rt+betatl*r(tg)))
      return
    end if
! 
  end function carbon
!
end module mod_co2
