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

! Unit test of Share/mod_capecin.F90 on two soundings: the analytic sounding
! of Weisman and Klemp (1982) with 14 g/kg of surface mixing ratio, which
! has a CAPE of about 2000-2500 J/kg and a negative lifted index, and a dry
! isothermal atmosphere at 250 K, without CAPE, whose lifted index is the
! difference from a dry adiabat, 250 - 250*(500/1000)**(R/cp) = 45 K.

program test_capecin
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil
  use mod_capecin
  use share_test_utils
  implicit none
  integer(ik4), parameter :: nk = 81
  real(rkx), pointer, contiguous, dimension(:,:,:) :: p => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: t => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: rh => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: qv => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: td => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: pi => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: q => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: th => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: thv => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: z => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: pr => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: tr => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: qr => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: cape => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: cin => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: li => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: t500 => null()
  integer(ik4) :: n

  call memory_init
  call getmem(p,1,1,1,2,1,nk,'p')
  call getmem(t,1,1,1,2,1,nk,'t')
  call getmem(rh,1,1,1,2,1,nk,'rh')
  call getmem(qv,1,1,1,2,1,nk,'qv')
  call getmem(td,1,1,1,2,1,nk,'td')
  call getmem(pi,1,1,1,2,1,nk,'pi')
  call getmem(q,1,1,1,2,1,nk,'q')
  call getmem(th,1,1,1,2,1,nk,'th')
  call getmem(thv,1,1,1,2,1,nk,'thv')
  call getmem(z,1,1,1,2,1,nk,'z')
  call getmem(pr,1,1,1,2,1,nk,'pr')
  call getmem(tr,1,1,1,2,1,nk,'tr')
  call getmem(qr,1,1,1,2,1,nk,'qr')
  call getmem(cape,1,1,1,2,'cape')
  call getmem(cin,1,1,1,2,'cin')
  call getmem(li,1,1,1,2,'li')
  call getmem(t500,1,1,1,2,'t500')

  call weisman_klemp(1)
  call isothermal(2)
  call getcape_new(1,1,1,2,nk,p,t,rh,cape,cin,td,pi,q,th,thv,z)
  call check(cape(1,1) > 1800.0_rkx .and. cape(1,1) < 2800.0_rkx, &
    'CAPE of the Weisman-Klemp sounding')
  call check(cin(1,1) >= 0.0_rkx .and. cin(1,1) < 100.0_rkx, &
    'CIN of the Weisman-Klemp sounding')
  call check(cape(1,2) == 0.0_rkx,'no CAPE in a dry isothermal atmosphere')

  ! otlift: levels from the top, the parcel from the lowest level.
  do n = 1, nk
    pr(1,:,n) = p(1,:,nk-n+1)
    tr(1,:,n) = t(1,:,nk-n+1)
    qr(1,:,n) = qv(1,:,nk-n+1)
  end do
  do n = 1, 2
    t500(1,n) = interp500(n)
  end do
  call otlift(li,tr,qr,pr,t500,1,1,1,2,nk)
  call check(li(1,1) < -3.0_rkx,'lifted index of the Weisman-Klemp sounding')
  call check(abs(li(1,2)-(250.0_rkx-250.0_rkx*0.5_rkx**rovcp)) &
    < 0.5_rkx,'lifted index of a dry isothermal atmosphere')

  call memory_destroy
  call test_end('test_capecin')

  contains

  ! Levels every 250 m from the surface (1000 hPa) to 20 km.
  subroutine weisman_klemp(n)
    implicit none
    integer(ik4), intent(in) :: n
    real(rkx), parameter :: th0 = 300.0_rkx, thtr = 343.0_rkx
    real(rkx), parameter :: ttr = 213.0_rkx, ztr = 12000.0_rkx
    real(rkx) :: zz, thz, rhz, pik, pikm1, tk, qs
    integer(ik4) :: k
    pikm1 = 1.0_rkx
    do k = 1, nk
      zz = 250.0_rkx*real(k-1,rkx)
      if ( zz <= ztr ) then
        thz = th0 + (thtr-th0)*(zz/ztr)**1.25_rkx
        rhz = 1.0_rkx - 0.75_rkx*(zz/ztr)**1.25_rkx
      else
        thz = thtr*exp(egrav*(zz-ztr)/(cpd*ttr))
        rhz = 0.25_rkx
      end if
      if ( k == 1 ) then
        pik = 1.0_rkx
      else
        pik = pikm1 - egrav*250.0_rkx/(cpd*0.5_rkx*(thz+th(1,n,k-1)))
      end if
      th(1,n,k) = thz
      p(1,n,k) = 100000.0_rkx*pik**(1.0_rkx/rovcp)
      tk = thz*pik
      t(1,n,k) = tk
      qs = ep2*es(tk)/(p(1,n,k)-es(tk))
      qv(1,n,k) = min(rhz*qs,0.014_rkx)
      rh(1,n,k) = qv(1,n,k)/qs
      pikm1 = pik
    end do
  end subroutine weisman_klemp

  subroutine isothermal(n)
    implicit none
    integer(ik4), intent(in) :: n
    integer(ik4) :: k
    do k = 1, nk
      p(1,n,k) = 100000.0_rkx*exp(-egrav*250.0_rkx*real(k-1,rkx) / &
                                  (rgas*250.0_rkx))
      t(1,n,k) = 250.0_rkx
      rh(1,n,k) = 0.001_rkx
      qv(1,n,k) = 0.001_rkx*ep2*es(250.0_rkx)/(p(1,n,k)-es(250.0_rkx))
    end do
  end subroutine isothermal

  ! Saturation vapour pressure (Pa), Bolton (1980)
  real(rkx) function es(tk)
    implicit none
    real(rkx), intent(in) :: tk
    es = 611.2_rkx*exp(17.67_rkx*(tk-273.15_rkx)/(tk-29.65_rkx))
  end function es

  ! Virtual temperature at 500 hPa, linear in ln p
  real(rkx) function interp500(n)
    implicit none
    integer(ik4), intent(in) :: n
    integer(ik4) :: k
    real(rkx) :: w, tv1, tv2
    interp500 = t(1,n,nk)
    do k = 1, nk-1
      if ( p(1,n,k) >= 50000.0_rkx .and. p(1,n,k+1) < 50000.0_rkx ) then
        w = log(p(1,n,k)/50000.0_rkx)/log(p(1,n,k)/p(1,n,k+1))
        tv1 = t(1,n,k)*(1.0_rkx+ep1*qv(1,n,k))
        tv2 = t(1,n,k+1)*(1.0_rkx+ep1*qv(1,n,k+1))
        interp500 = (1.0_rkx-w)*tv1 + w*tv2
      end if
    end do
  end function interp500

end program test_capecin
