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

module mod_bats_co2
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_bats_param
  use mod_bats_leaftemp
  use mod_bats_internal

  implicit none

  private

  public :: co2

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
!            cari  = rate of carbon uptake by plants
!            apbm  = accumulated primary biomass
!              ra  = leaf aerodynamic resistance factor
!             rap  = leaf boundary layer resistance to co2
!           resps  = soil respiration
!             rmp  = physical mesophyllic resistance
!              rs  = stomatal resistance
!             rsp  = stomatal resistance to co2
!              rt  = total mechanical resistance to co2
!      resp(i) = carbon in soil and in veg (kg c / m**2 / sec)
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
    implicit none
    integer(ik4) :: i
    real(rkx) :: rap , resps , rsp , rt , rcar , cari , apbm
    real(rkx) , parameter :: rmp = 800.0_rkx

    apbm = d_zero
    do i = ilndbeg , ilndend
      if ( sigf(i) > minsigf ) then
        rsp = lftrs(i)*1.7_rkx
        rap = lftra(i)*1.5_rkx
        rt = rsp + rap + rmp
        rcar = carbon(swsi(i)*rlai(i),tlef(i),rt,tgrd(i),xlai(i),xlsai(i))
        cari = sigf(i)*xlsai(i)*fdry(i)*rcar
        apbm = apbm + cari*dtbat
        if ( apbm < d_zero ) apbm = d_zero
        resps = 0.7e-7_rkx*resp(i)*exp(0.1_rkx*(tgrd(i)-300.0_rkx)) * &
                min(d_one,ssw(i)/(0.6_rkx*gwmx0(i)))
        resp(i) = resp(i) + (cari-resps)*dtbat
        if ( resp(i) < d_zero ) resp(i) = d_zero
      end if
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
  pure real(rkx) function carbon(vf,t,rm,tg,xlai,xlsai)
    implicit none
    real(rkx) , intent(in) :: rm , t , tg , vf , xlai , xlsai
    real(rkx) :: ab , ac , al , alphtl , b , bc , betatl , cco2 ,   &
               cco2i , ccold , gt , p , pm , pml , rt , w , &
               wd , wp , xk , xkb , xl
    integer(ik4) :: it
    !
    !====================================================================
    ! this fn is not in si - vf (radn) and rm (res) handed in are in si
    ! result of fn (carbon) is handed back in si (kg/m**2/s)
    !======================================================================
    !
    p = d_zero
    rt = r(t)
    ! alphtl=sai/lai     betatl=rai/lai  (rai=root area index)
    alphtl = (xlsai-xlai)/xlai
    betatl = d_half
    if ( vf < d_two ) then
      ! nighttime maintenance respiration only
      carbon = -0.36e-8_rkx*((d_one+alphtl)*0.877_rkx*rt+ &
                 betatl*(r(tg)-0.123_rkx*rt))
    else
      ! convert lambda less than 0.7 micron solar into photon units
      ! light intensity (e/m2)
      xl = 4.6e-3_rkx*vf
      ! co2 external concentration(mm/m3)
      ! 335ppm/v
      cco2 = 13.5_rkx
      ! initial guess for co2 concentration inside chloroplast
      cco2i = cco2
      ! co2 half max in absence of oxygen  (mm/m3)
      xk = d_half
      ! oxygen inhibition factor
      b = 3.56_rkx
      xkb = xk*b
      ! maximum temperature optimum light saturated photosynthesis
      pml = 0.050_rkx
      ! quantum yield(mm/e)
      al = 0.05_rkx
      gt = g(t,320.0_rkx,4.0e3_rkx)
      ! maximum photosynthesis
      pm = e(xl,al,pml*gt)
      ! iterate
      do it = 1 , 30
        ! photorespiration
        wp = pm/(d_one+0.4_rkx*(d_one+cco2i/xk))
        ! total respiration
        ! dark respiration within daytime leaves
        wd = 3.0e-4_rkx*rt + 0.14_rkx*p
        w = wp + wd
        ! carbon uptake factors
        ac = cco2 + xkb + rm*(pm-w)
        bc = d_four*rm*(cco2*(pm-w)-xkb*w)
        ab = sqrt(ac**2-bc)
        p = d_half*(ac-ab)/rm
        ccold = cco2i
        cco2i = cco2 - p*rm
        if ( abs(cco2i-ccold) <= 0.05_rkx .and. it > 9 ) exit
      end do
      ! respiration outside leaves
      carbon = 1.2e-5_rkx*((d_one-0.14_rkx*(alphtl+betatl)) * &
               p-3.0e-4_rkx*(alphtl*rt+betatl*r(tg)))
    end if

    contains

      pure real(rkx) function g(t,tmx,sl)
        implicit none
        real(rkx) , intent(in) :: t , tmx , sl
        g = exp(sl*(d_one/tmx-d_one/t)) / &
            (d_one+(exp(sl*(d_one/tmx-d_one/t)*6.0_rkx)))*5.0e-3_rkx*t
      end function g
      ! temperature dependence of dark respiration
      pure real(rkx) function r(t)
        implicit none
        real(rkx) , intent(in) :: t
        r = exp(30.0_rkx-9.0e3_rkx/t)
      end function r
      ! light dependence of photosynthesis
      pure real(rkx) function e(xl,a,pml)
        implicit none
        real(rkx) , intent(in) :: xl , a , pml
        e = a*xl/sqrt(d_one+(a*xl/pml)**2)
      end function e

  end function carbon
!
end module mod_bats_co2
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
