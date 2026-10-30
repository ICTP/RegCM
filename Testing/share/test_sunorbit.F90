!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Unit test of Share/mod_sunorbit.F90: the orbital parameters for 1950
! against the reference values of Berger (1978), the
! declination at the solstices and equinoxes, the earth-sun distance factor
! at perihelion and aphelion and its annual mean, and the cosine of the
! zenith angle at noon and midnight.

program test_sunorbit
  use mod_intkinds
  use mod_realkinds
  use mod_constants, only : mathpi, degrad
  use mod_dynparam, only : dayspy, vernal_equinox
  use mod_sunorbit
  use share_test_utils
  implicit none
  real(rk8) :: eccen, obliq, mvelp, obliqr, lambm0, mvelpp
  real(rk8) :: delta, eccf, day, emax, emin, dmax, dmin, esum, d
  integer(ik4) :: n, dmaxday, emaxday, eminday

  dayspy = 365.0_rkx
  vernal_equinox = 80.5_rkx

  call orb_params(1950,eccen,obliq,mvelp,obliqr,lambm0,mvelpp)
  call check(abs(eccen-0.016724_rk8) < 1.0e-6_rk8,'eccentricity in 1950')
  call check(abs(obliq-23.446_rk8) < 1.0e-3_rk8,'obliquity in 1950')
  call check(abs(mvelp-102.04_rk8) < 0.01_rk8, &
    'longitude of perihelion in 1950')
  call check(abs(obliqr-obliq*degrad) < 1.0e-14_rk8, &
    'obliquity in radians')

  ! Along the year: extremes of declination and distance factor.
  dmax = -1.0_rk8
  dmin = 1.0_rk8
  emax = 0.0_rk8
  emin = 2.0_rk8
  esum = 0.0_rk8
  dmaxday = 0
  emaxday = 0
  eminday = 0
  do n = 0, 3649
    day = 1.0_rk8 + real(n,rk8)*0.1_rk8
    call orb_decl(day,eccen,mvelpp,lambm0,obliqr,delta,eccf)
    if ( delta > dmax ) then
      dmax = delta
      dmaxday = int(day)
    end if
    dmin = min(dmin,delta)
    if ( eccf > emax ) then
      emax = eccf
      emaxday = int(day)
    end if
    if ( eccf < emin ) then
      emin = eccf
      eminday = int(day)
    end if
    esum = esum + eccf
  end do
  call check(abs(dmax-obliqr) < 1.0e-4_rk8 .and. &
             abs(dmin+obliqr) < 1.0e-4_rk8, &
    'declination between minus and plus the obliquity')
  call check(abs(dmaxday-172) <= 2, &
    'maximum declination at the June solstice')
  call check(abs(emax-((1.0_rk8+eccen)/(1.0_rk8-eccen**2))**2) &
    < 1.0e-6_rk8 .and. &
    abs(emin-((1.0_rk8-eccen)/(1.0_rk8-eccen**2))**2) < 1.0e-6_rk8, &
    'distance factor at perihelion and aphelion')
  call check(emaxday >= 1 .and. emaxday <= 6 .and. &
             abs(eminday-185) <= 3,'perihelion early in January')
  call check(abs(esum/3650.0_rk8-1.0_rk8/sqrt(1.0_rk8-eccen**2)) &
    < 1.0e-4_rk8, &
    'annual mean of the distance factor')
  call orb_decl(real(vernal_equinox,rk8),eccen,mvelpp,lambm0,obliqr, &
                delta,eccf)
  call check(abs(delta) < 2.0e-3_rk8,'zero declination at the equinox')

  ! Zenith angle: day 80.5 is noon at Greenwich at the equinox.
  d = 0.0_rk8
  call check(abs(orb_cosz(80.5_rk8,0.0_rk8,0.0_rk8,d)-1.0_rk8) &
    < 1.0e-12_rk8, &
    'sun at the zenith at noon on the equator')
  call check(abs(orb_cosz(80.0_rk8,0.0_rk8,0.0_rk8,d)+1.0_rk8) &
    < 1.0e-12_rk8, &
    'sun at the nadir at midnight')
  call check(abs(orb_cosz(80.25_rk8,0.0_rk8,90.0_rk8,d)-1.0_rk8) &
    < 1.0e-12_rk8,'noon 6 hours earlier at 90E')
  call check(abs(orb_cosz(80.5_rk8,60.0_rk8,0.0_rk8,d)- &
             cos(60.0_rk8*degrad)) < 1.0e-12_rk8,'noon at 60N')
  call check(abs(real(orb_cosz(80.5_rk4,60.0_rk4,0.0_rk4,d),rk8)- &
             0.5_rk8) < 1.0e-6_rk8,'single precision')

  call test_end('test_sunorbit')
end program test_sunorbit
