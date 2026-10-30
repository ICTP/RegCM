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

! Unit test of Share/mod_vertint.F90 on exact profiles: interpolation in
! ln p of a field linear in ln p (intlog, intv2, intv3) with the standard
! lapse rate extrapolation below the lowest level, linear interpolation in
! p (intlinreg), the layer mean temperature of intgtb, and the surface
! pressure of a standard atmosphere at any surface height (intzps).

program test_vertint
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_vertint
  use mod_memutil
  use share_test_utils
  implicit none
  integer(ik4), parameter :: km = 8, kp = 5, nk = 12, nz = 7, ng = 6
  real(rk8), parameter :: p0 = 101325.0_rk8
  real(rk8), dimension(km) :: pl, sig
  real(rk8), dimension(1,1,km) :: f3, p3
  real(rk8), dimension(1,1,kp) :: fp
  real(rk8), dimension(kp) :: pt, ex
  real(rk8), dimension(1,1) :: ps8
  real(rk8) :: ptop
  real(rkx), dimension(km) :: sccm
  real(rkx), dimension(kp) :: srcm, exr
  real(rkx), dimension(1,1,km) :: fc
  real(rkx), dimension(1,1,kp) :: fr
  real(rkx), dimension(1,1) :: psr, fs
  real(rkx) :: pss, pst
  real(rk8), pointer, contiguous, dimension(:,:,:) :: lf => null()
  real(rk8), pointer, contiguous, dimension(:,:,:) :: lfp => null()
  real(rk8), pointer, contiguous, dimension(:,:,:) :: lp3 => null()
  real(rk8), pointer, contiguous, dimension(:,:) :: lps => null()
  real(rk8), pointer, contiguous, dimension(:) :: lpl => null()
  real(rkx), dimension(ng) :: sg
  real(rkx), dimension(1,1,ng) :: tg, zg
  real(rkx), dimension(1,1) :: zr, tl, pa, za
  real(rkx), dimension(nk), parameter :: plev = &
    [1000.0_rkx, 975.0_rkx, 950.0_rkx, 925.0_rkx, 900.0_rkx, 875.0_rkx, &
      850.0_rkx, 800.0_rkx, 700.0_rkx, 600.0_rkx, 500.0_rkx, 300.0_rkx]
  real(rkx), dimension(nz), parameter :: zs = &
    [-50.0_rkx, 100.0_rkx, 300.0_rkx, 800.0_rkx, 1500.0_rkx, &
     2500.0_rkx, 4000.0_rkx]
  real(rkx), dimension(3), parameter :: t0s = &
    [288.15_rkx, 303.15_rkx, 265.15_rkx]
  real(rkx), dimension(3), parameter :: zpstol = &
    [0.02_rkx, 0.5_rkx, 0.02_rkx]
  real(rkx), pointer, contiguous, dimension(:,:) :: ps => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: zrp => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: lat => null()
  real(rkx), pointer, contiguous, dimension(:) :: sz => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: tz => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: zz => null()
  real(rk8) :: t0, err
  integer(ik4) :: k, n, m
  character(len=80) :: msg

  call memory_init

  ! A field linear in ln p on levels from 1000 to 100 hPa. Below the
  ! lowest level the extrapolation follows the standard lapse rate.
  pl = [1000.0_rk8, 925.0_rk8, 850.0_rk8, 700.0_rk8, 500.0_rk8, &
         300.0_rk8, 200.0_rk8, 100.0_rk8]
  f3(1,1,:) = 250.0_rk8 + 20.0_rk8*log(pl)
  p3(1,1,:) = pl
  pt = [960.0_rk8, 780.0_rk8, 400.0_rk8, 150.0_rk8, 1030.0_rk8]
  ex(1:4) = 250.0_rk8 + 20.0_rk8*log(pt(1:4))
  ex(5) = 0.5_rk8*(f3(1,1,1)+f3(1,1,2)) * &
          exp(rgas*lrate/egrav*log(pt(5)/pl(1)))
  call intlog(fp,f3,p3,1,1,km,pt,kp)
  call check(maxval(abs(fp(1,1,:)-ex)) < 1.0e-9_rk8, &
    'intlog on pressure levels')
  ptop = 100.0_rk8
  ps8 = 1000.0_rk8
  sig = (pl-ptop)/(ps8(1,1)-ptop)
  call intlog(fp,f3,ps8,sig,ptop,1,1,km,pt,kp)
  call check(maxval(abs(fp(1,1,:)-ex)) < 1.0e-9_rk8, &
    'intlog on sigma levels')

  ! intv2 and intv3: sigma input (pss, pst in kPa) to model sigma levels
  ! with p* = 95 kPa and ptop = 5 kPa.
  pss = real((pl(1)-pl(km))/10.0_rk8,rkx)
  pst = real(pl(km)/10.0_rk8,rkx)
  sccm = real((pl-pl(km))/(pl(1)-pl(km)),rkx)
  fc(1,1,:) = real(f3(1,1,:),rkx)
  psr = 95.0_rkx
  do k = 1, kp
    srcm(k) = real((pt(k)/10.0_rk8-5.0_rk8)/95.0_rk8,rkx)
  end do
  call intv2(fr,fc,psr,srcm,pss,sccm,5.0_rkx,pst,1,1,kp,km)
  exr = real(ex,rkx)
  exr(5) = real(f3(1,1,1)*(-egrav/lrate-rgas*0.5_rk8*log(pt(5)/pl(1))) / &
                (-egrav/lrate+rgas*0.5_rk8*log(pt(5)/pl(1))),rkx)
  call check(maxval(abs(fr(1,1,:)-exr)) < 1.0e-4_rkx,'intv2')
  call intv3(fs,fc,psr,pss,sccm,5.0_rkx,pst,1,1,km)
  call check(abs(fs(1,1)-real(250.0_rk8+20.0_rk8*log(1000.0_rk8),rkx)) &
    < 1.0e-4_rkx,'intv3 at the surface')

  ! intlinreg: a field linear in p, levels ordered both ways. Outside the
  ! levels the nearest value is used.
  call getmem(lf,1,1,1,1,1,km,'lf')
  call getmem(lfp,1,1,1,1,1,kp,'lfp')
  call getmem(lp3,1,1,1,1,1,kp,'lp3')
  call getmem(lps,1,1,1,1,'lps')
  call getmem(lpl,1,km,'lpl')
  lps = 1000.0_rk8
  lp3(1,1,:) = [960.0_rk8, 780.0_rk8, 400.0_rk8, 150.0_rk8, 50.0_rk8]
  do n = 1, 2
    if ( n == 1 ) then
      lpl = pl
    else
      lpl = pl(km:1:-1)
    end if
    lf(1,1,:) = 3.0_rk8*lpl - 100.0_rk8
    call intlinreg(lfp,lf,lps,lpl,1,1,1,1,km,lp3,kp)
    err = maxval(abs(lfp(1,1,1:4)-(3.0_rk8*lp3(1,1,1:4)-100.0_rk8)))
    err = max(err,abs(lfp(1,1,5)-(3.0_rk8*100.0_rk8-100.0_rk8)))
    write(msg,'(a,i0)') 'intlinreg, level order ', n
    call check(err < 1.0e-9_rk8,trim(msg))
  end do

  ! intgtb: T linear in z, the layer mean is the mean of the two ends.
  zg(1,1,:) = [100.0_rkx, 400.0_rkx, 800.0_rkx, 1500.0_rkx, &
               3000.0_rkx, 6000.0_rkx]
  tg(1,1,:) = 290.0_rkx - 0.0065_rkx*zg(1,1,:)
  sg = [1.0_rkx, 0.95_rkx, 0.9_rkx, 0.8_rkx, 0.6_rkx, 0.3_rkx]
  zr = 500.0_rkx
  call intgtb(pa,za,tl,zr,tg,zg,100.0_rkx,sg,0.0_rkx,1,1,ng)
  call check(abs(tl(1,1)-0.5_rkx*((290.0_rkx-0.0065_rkx*500.0_rkx) + &
             (290.0_rkx-0.0065_rkx*800.0_rkx))) < 1.0e-4_rkx, &
    'intgtb layer mean temperature')

  ! intzps: standard atmospheres with three sea level temperatures. With
  ! a warm atmosphere the ECMWF limit on T* changes the extrapolation
  ! below the lowest level.
  call getmem(ps,1,1,1,1,'ps')
  call getmem(zrp,1,1,1,1,'zrp')
  call getmem(lat,1,1,1,1,'lat')
  call getmem(sz,1,nk,'sz')
  call getmem(tz,1,1,1,1,1,nk,'tz')
  call getmem(zz,1,1,1,1,1,nk,'zz')
  lat = 0.0_rkx
  pss = (plev(1)-plev(nk))/10.0_rkx
  pst = plev(nk)/10.0_rkx
  sz = (plev-plev(nk))/(plev(1)-plev(nk))
  do m = 1, 3
    t0 = real(t0s(m),rk8)
    do k = 1, nk
      zz(1,1,k) = real((t0/lrate) * &
        (1.0_rk8-(real(plev(k),rk8)*100.0_rk8/p0)**(rgas*lrate/egrav)),rkx)
      tz(1,1,k) = real(t0-lrate*real(zz(1,1,k),rk8),rkx)
    end do
    err = 0.0_rk8
    do n = 1, nz
      zrp = zs(n)
      call intzps(ps,zrp,tz,zz,pss,sz,pst,lat,1.0_rkx,365.0_rkx,1,1,nk)
      err = max(err,abs(real(ps(1,1),rk8)*10.0_rk8 - p0 * &
        ((t0-lrate*real(zs(n),rk8))/t0)**(egrav/(rgas*lrate))/100.0_rk8))
    end do
    write(msg,'(a,f6.2,a)') 'intzps, sea level temperature ', t0, ' K'
    call check(err < real(zpstol(m),rk8),trim(msg))
  end do

  call memory_destroy
  call test_end('test_vertint')
end program test_vertint
