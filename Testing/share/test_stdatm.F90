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

! Unit test of Share/mod_stdatm.F90 through stdatm_val: values at the
! table levels, the zones at the equator and at the poles, continuity at
! 45 degrees, the seasons in the two hemispheres, the air density against
! p/(RT) and the decrease of water vapour with height.

program test_stdatm
  use mod_intkinds
  use mod_realkinds
  use mod_stdatm
  use share_test_utils
  implicit none
  real(rk8), parameter :: rd = 287.05_rk8
  real(rk8), parameter :: dpy = 365.0_rk8
  real(rk8), dimension(5), parameter :: lats = &
    [0.0_rk8, 30.0_rk8, 50.0_rk8, 70.0_rk8, 90.0_rk8]
  real(rk8) :: p, t, rho, q, qold, err, lat, day, v1, v2
  integer(ik4) :: k, n, m, nbad

  ! Equator: the tropical profile. North pole: mean of the subarctic
  ! winter and summer profiles.
  call check(abs(stdatm_val(0.0_rk8,1013.0_rk8,istdatm_tempk)- &
                 300.0_rk8) < 1.0e-9_rk8,'tropical surface temperature')
  call check(abs(stdatm_val(0.0_rk8,286.0_rk8,istdatm_tempk)- &
                 237.0_rk8) < 1.0e-9_rk8,'tropical 10 km temperature')
  call check(abs(stdatm_val(90.0_rk8,2000.0_rk8,istdatm_tempk)- &
                 0.5_rk8*(257.1_rk8+287.0_rk8)) < 1.0e-9_rk8, &
    'polar surface temperature, mean of the seasons')

  ! Continuity across 45 degrees.
  err = 0.0_rk8
  do k = 1, 40
    p = 1000.0_rk8*exp(-0.1_rk8*real(k,rk8))
    err = max(err,abs(stdatm_val(44.9999_rk8,p,istdatm_tempk)- &
                      stdatm_val(45.0001_rk8,p,istdatm_tempk)))
  end do
  call check(err < 1.0e-2_rk8,'continuous at 45 degrees')

  ! The two hemispheres have opposite seasons.
  err = 0.0_rk8
  do n = 2, 4
    do m = 0, 11
      day = 1.0_rk8 + real(m,rk8)*dpy/12.0_rk8
      v1 = stdatm_val(day,dpy,lats(n),500.0_rk8,istdatm_tempk)
      v2 = stdatm_val(day+dpy*0.5_rk8,dpy,-lats(n),500.0_rk8, &
                      istdatm_tempk)
      err = max(err,abs(v1-v2))
    end do
  end do
  call check(err < 1.0e-6_rk8,'opposite seasons in the two hemispheres')
  v1 = stdatm_val(15.0_rk8,dpy,60.0_rk8,1000.0_rk8,istdatm_tempk)
  v2 = stdatm_val(196.0_rk8,dpy,60.0_rk8,1000.0_rk8,istdatm_tempk)
  call check(v1 < v2,'January colder than July at 60N')
  v1 = stdatm_val(15.0_rk8,dpy,-60.0_rk8,1000.0_rk8,istdatm_tempk)
  v2 = stdatm_val(196.0_rk8,dpy,-60.0_rk8,1000.0_rk8,istdatm_tempk)
  call check(v1 > v2,'January warmer than July at 60S')

  ! Air density consistent with pressure and temperature at any latitude,
  ! season and height up to 3 hPa. Between levels 5 km apart the
  ! interpolation, linear in ln p, overestimates the density by up to 9 %;
  ! the errors of the table that this check found were 30 %.
  nbad = 0
  do n = 1, 5
    do m = 0, 3
      day = 1.0_rk8 + real(m,rk8)*dpy/4.0_rk8
      do k = 1, 49
        p = 1013.0_rk8*exp(-0.12_rk8*real(k-1,rk8))
        t = stdatm_val(day,dpy,lats(n),p,istdatm_tempk)
        rho = stdatm_val(day,dpy,lats(n),p,istdatm_airdn)
        if ( abs(rho/(p*1.0e5_rk8/(rd*t))-1.0_rk8) > 0.10_rk8 ) &
          nbad = nbad + 1
      end do
    end do
  end do
  call check(nbad == 0,'density consistent with p/(RT)')

  ! Water vapour decreases with height in the troposphere.
  nbad = 0
  do n = 1, 5
    do m = 0, 3
      day = 1.0_rk8 + real(m,rk8)*dpy/4.0_rk8
      lat = lats(n)
      qold = 1.0e30_rk8
      do k = 1, 15
        p = 1000.0_rk8 - 50.0_rk8*real(k-1,rk8)
        q = stdatm_val(day,dpy,lat,p,istdatm_qdens)
        if ( q > qold*1.0001_rk8 ) nbad = nbad + 1
        qold = q
      end do
    end do
  end do
  call check(nbad == 0,'water vapour decreases up to 300 hPa')

  ! Above the top of the table: the value at the top.
  call check(abs(stdatm_val(0.0_rk8,0.1_rk8,istdatm_tempk)-270.0_rk8) &
    < 1.0e-9_rk8,'constant above the table')

  call test_end('test_stdatm')
end program test_stdatm
