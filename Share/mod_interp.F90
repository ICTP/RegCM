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

module mod_interp

  !
  ! One dimensional interpolation with a tension spline (interp1d), used
  ! for the vertical interpolation (mod_vertint) and for the aerosol
  ! profiles.
  !
  use mod_intkinds
  use mod_realkinds

  implicit none

  private

  public :: interp1d

  interface interp1d
    module procedure interp1d_r4
    module procedure interp1d_r8
    module procedure interp1d_r4_new
    module procedure interp1d_r8_new
    module procedure single_interp1d_r4
    module procedure single_interp1d_r8
  end interface interp1d

  interface merge_coincident
    module procedure merge_coincident_k4
    module procedure merge_coincident_k8
  end interface merge_coincident

  contains

  pure subroutine single_interp1d_r4(xi,g,xo,f,alfa,ex1,ex2)
    !$acc routine seq
    implicit none
    real(rk4), dimension(:), intent(in) :: xi, g
    real(rk4), intent(in) :: xo
    real(rk4), intent(out) :: f
    real(rk4), intent(in) :: alfa
    real(rk4), intent(in) :: ex1
    real(rk4), intent(in) :: ex2

    real(rk4) :: fmm, fpp, xmm, xpp
    real(rk4) :: fm, xm, fp, xp, delx, delxp, delxm
    real(rk4) :: delx1, delx2, delxs, delx1s, delx2s
    real(rk4) :: spl, clin

    real(rk4), dimension(size(xi)) :: zi, zg
    integer(ik4) :: npi, k, jj, ir

    npi = size(xi)
    if ( xi(1) >= xi(npi) ) then
      do k = 1, npi
        zi(k) = xi(npi-k+1)
        zg(k) = g(npi-k+1)
      end do
    else
      zi(:) = xi(:)
      zg(:) = g(:)
    end if

    ! Merge coincident or non increasing coordinates
    call merge_coincident(zi,zg,npi)
    if ( npi == 1 ) then
      f = zg(1)
      return
    end if

    if ( xo < zi(1) ) then
      f = zg(1) + ex1*(zg(1)-zg(2))/(zi(1)-zi(2)) * (xo-zi(1))
      return
    else if ( xo >= zi(npi) ) then
      f = zg(npi) + ex2*(zg(npi)-zg(npi-1))/(zi(npi)-zi(npi-1)) * &
           (xo-zi(npi))
      return
    end if
    if ( npi == 2 ) then
      f = (zg(1)*(zi(2)-xo) + zg(2)*(xo-zi(1)))/(zi(2)-zi(1))
      return
    end if
    ir = 0
    do jj = 1, npi
      if ( xo >= zi(jj) ) ir = ir + 1
    end do
    if ( ir == 1 ) then
      fmm = 2.0_rk4 * zg(1) - zg(2)
      xmm = 2.0_rk4 * zi(1) - zi(2)
      fpp = zg(ir+2)
      xpp = zi(ir+2)
    else if ( ir == (npi-1) ) then
      fpp = 2.0_rk4 * zg(npi) - zg(npi-1)
      xpp = 2.0_rk4 * zi(npi) - zi(npi-1)
      fmm = zg(ir-1)
      xmm = zi(ir-1)
    else
      fmm = zg(ir-1)
      xmm = zi(ir-1)
      fpp = zg(ir+2)
      xpp = zi(ir+2)
    end if
    fm     = zg(ir)
    xm     = zi(ir)
    fp     = zg(ir+1)
    xp     = zi(ir+1)
    delx   = xp - xm
    delxp  = xpp - xp
    delxm  = xm - xmm
    delx1  = xo - xm
    delx2  = xp - xo
    delxs  = delx**2
    delx1s = delx1**2
    delx2s = delx2**2
    spl = fm*(delx2/delx + delx1*delx2s/(delxs*delxm) - delx1s*     &
          delx2/((delx+delxp)*delxs)) + fp*(delx1/delx +            &
          delx1s*delx2/(delxs*delxp) - delx1*delx2s/((delx+delxm)*  &
          delxs)) - fmm * delx1*delx2s/((delx+delxm)*delx*delxm) -  &
          fpp * delx1s*delx2/((delx+delxp)*delx*delxp)
    clin = (fm*delx2 + fp*delx1)/delx
    f = alfa*clin + (1.0_rk4-alfa)*spl
  end subroutine single_interp1d_r4

  pure subroutine interp1d_r4(xi,g,xo,f,alfa,ex1,ex2)
    !$acc routine seq
    implicit none
    !  Input:  function g defined at irregular but strictly monotonic xi
    !  Output: f interpolated values at arbitrary coordinates xo
    real(rk4), dimension(:), intent(in) :: xi, xo, g
    real(rk4), dimension(:), intent(out) :: f
    ! alfa: spline tension parameter, comprised between 0 and 1:
    ! if alfa=1, pure linear interpolation; if alfa=0, pure spline
    real(rk4), intent(in) :: alfa
    ! ex1: param. determining extrapolation for x < xi(1)
    ! ex2: param. determining extrapolation for x > xi(npi)
    ! if ex1=0 or ex2=0, constant value extrapolation is used at extreme
    ! if ex1=1 or ex2=1, linear extrapolation is used at corresponding extreme
    ! intermediate values of ex1 and ex2 give intermediate extrapolation values
    real(rk4), intent(in) :: ex1
    real(rk4), intent(in) :: ex2

    real(rk4) :: fmm, fpp, xmm, xpp
    real(rk4) :: fm, xm, fp, xp, delx, delxp, delxm
    real(rk4) :: delx1, delx2, delxs, delx1s, delx2s
    real(rk4) :: spl, clin

    real(rk4), dimension(size(xi)) :: zi, zg
    integer(ik4) :: npi, npo
    integer(ik4) :: k, j, jj, ir

    npi = size(xi)
    npo = size(xo)
    if ( xi(1) >= xi(npi) ) then
      do k = 1, npi
        zi(k) = xi(npi-k+1)
        zg(k) = g(npi-k+1)
      end do
    else
      zi(:) = xi(:)
      zg(:) = g(:)
    end if

    ! Merge coincident or non increasing coordinates
    call merge_coincident(zi,zg,npi)
    if ( npi == 1 ) then
      f(:) = zg(1)
      return
    end if

    do j = 1, npo
      ! 2 cases of extrapolation
      if ( xo(j) < zi(1) ) then
        f(j) = zg(1) + ex1*(zg(1)-zg(2))/(zi(1)-zi(2)) * (xo(j)-zi(1))
        cycle
      else if ( xo(j) >= zi(npi) ) then
        f(j) = zg(npi) + ex2*(zg(npi)-zg(npi-1))/(zi(npi)-zi(npi-1)) * &
               (xo(j)-zi(npi))
        cycle
      end if
      if ( npi == 2 ) then
        f(j) = (zg(1)*(zi(2)-xo(j)) + zg(2)*(xo(j)-zi(1))) / &
               (zi(2)-zi(1))
        cycle
      end if
      ir = 0
      ! ir is a reference index determining the interpolation interval
      ! The interpolation expression is applied also if xo = zi(j)
      do jj = 1, npi
        if ( xo(j) >= zi(jj) ) ir = ir + 1
      end do
      if ( ir == 1 ) then
        fmm = 2.0_rk4 * zg(1) - zg(2)
        xmm = 2.0_rk4 * zi(1) - zi(2)
        fpp = zg(ir+2)
        xpp = zi(ir+2)
      else if ( ir == (npi-1) ) then
        fpp = 2.0_rk4 * zg(npi) - zg(npi-1)
        xpp = 2.0_rk4 * zi(npi) - zi(npi-1)
        fmm = zg(ir-1)
        xmm = zi(ir-1)
      else
        fmm = zg(ir-1)
        xmm = zi(ir-1)
        fpp = zg(ir+2)
        xpp = zi(ir+2)
      end if
      fm     = zg(ir)
      xm     = zi(ir)
      fp     = zg(ir+1)
      xp     = zi(ir+1)
      delx   = xp - xm
      delxp  = xpp - xp
      delxm  = xm - xmm
      delx1  = xo(j) - xm
      delx2  = xp - xo(j)
      delxs  = delx**2
      delx1s = delx1**2
      delx2s = delx2**2
      !  Spline contribution to interpolation
      spl = fm*(delx2/delx + delx1*delx2s/(delxs*delxm) - delx1s*     &
            delx2/((delx+delxp)*delxs)) + fp*(delx1/delx +            &
            delx1s*delx2/(delxs*delxp) - delx1*delx2s/((delx+delxm)*  &
            delxs)) - fmm * delx1*delx2s/((delx+delxm)*delx*delxm) -  &
            fpp * delx1s*delx2/((delx+delxp)*delx*delxp)
      !  Linear interpolation contribution
      clin = (fm*delx2 + fp*delx1)/delx
      !  Final interpolation combined using alfa
      f(j) = alfa*clin + (1.0_rk4-alfa)*spl
    end do
  end subroutine interp1d_r4

  pure subroutine single_interp1d_r8(xi,g,xo,f,alfa,ex1,ex2)
    !$acc routine seq
    implicit none
    real(rk8), dimension(:), intent(in) :: xi, g
    real(rk8), intent(in) :: xo
    real(rk8), intent(out) :: f
    real(rk8), intent(in) :: alfa
    real(rk8), intent(in) :: ex1
    real(rk8), intent(in) :: ex2

    real(rk8) :: fmm, fpp, xmm, xpp
    real(rk8) :: fm, xm, fp, xp, delx, delxp, delxm
    real(rk8) :: delx1, delx2, delxs, delx1s, delx2s
    real(rk8) :: spl, clin

    real(rk8), dimension(size(xi)) :: zi, zg
    integer(ik4) :: npi, k, jj, ir

    npi = size(xi)
    if ( xi(1) >= xi(npi) ) then
      do k = 1, npi
        zi(k) = xi(npi-k+1)
        zg(k) = g(npi-k+1)
      end do
    else
      zi(:) = xi(:)
      zg(:) = g(:)
    end if

    ! Merge coincident or non increasing coordinates
    call merge_coincident(zi,zg,npi)
    if ( npi == 1 ) then
      f = zg(1)
      return
    end if

    if ( xo < zi(1) ) then
      f = zg(1) + ex1*(zg(1)-zg(2))/(zi(1)-zi(2)) * (xo-zi(1))
      return
    else if ( xo >= zi(npi) ) then
      f = zg(npi) + ex2*(zg(npi)-zg(npi-1))/(zi(npi)-zi(npi-1)) * &
           (xo-zi(npi))
      return
    end if
    if ( npi == 2 ) then
      f = (zg(1)*(zi(2)-xo) + zg(2)*(xo-zi(1)))/(zi(2)-zi(1))
      return
    end if
    ir = 0
    do jj = 1, npi
      if ( xo >= zi(jj) ) ir = ir + 1
    end do
    if ( ir == 1 ) then
      fmm = 2.0_rk8 * zg(1) - zg(2)
      xmm = 2.0_rk8 * zi(1) - zi(2)
      fpp = zg(ir+2)
      xpp = zi(ir+2)
    else if ( ir == (npi-1) ) then
      fpp = 2.0_rk8 * zg(npi) - zg(npi-1)
      xpp = 2.0_rk8 * zi(npi) - zi(npi-1)
      fmm = zg(ir-1)
      xmm = zi(ir-1)
    else
      fmm = zg(ir-1)
      xmm = zi(ir-1)
      fpp = zg(ir+2)
      xpp = zi(ir+2)
    end if
    fm     = zg(ir)
    xm     = zi(ir)
    fp     = zg(ir+1)
    xp     = zi(ir+1)
    delx   = xp - xm
    delxp  = xpp - xp
    delxm  = xm - xmm
    delx1  = xo - xm
    delx2  = xp - xo
    delxs  = delx**2
    delx1s = delx1**2
    delx2s = delx2**2
    spl = fm*(delx2/delx + delx1*delx2s/(delxs*delxm) - delx1s*     &
          delx2/((delx+delxp)*delxs)) + fp*(delx1/delx +            &
          delx1s*delx2/(delxs*delxp) - delx1*delx2s/((delx+delxm)*  &
          delxs)) - fmm * delx1*delx2s/((delx+delxm)*delx*delxm) -  &
          fpp * delx1s*delx2/((delx+delxp)*delx*delxp)
    clin = (fm*delx2 + fp*delx1)/delx
    f = alfa*clin + (1.0_rk8-alfa)*spl
  end subroutine single_interp1d_r8

  pure subroutine interp1d_r8(xi,g,xo,f,alfa,ex1,ex2)
    !$acc routine seq
    implicit none
    real(rk8), dimension(:), intent(in) :: xi, xo, g
    real(rk8), dimension(:), intent(out) :: f
    real(rk8), intent(in) :: alfa
    real(rk8), intent(in) :: ex1
    real(rk8), intent(in) :: ex2

    real(rk8) :: fmm, fpp, xmm, xpp
    real(rk8) :: fm, xm, fp, xp, delx, delxp, delxm
    real(rk8) :: delx1, delx2, delxs, delx1s, delx2s
    real(rk8) :: spl, clin

    real(rk8), dimension(size(xi)) :: zi, zg
    integer(ik4) :: npi, npo
    integer(ik4) :: k, j, jj, ir

    npi = size(xi)
    npo = size(xo)
    if ( xi(1) >= xi(npi) ) then
      do k = 1, npi
        zi(k) = xi(npi-k+1)
        zg(k) = g(npi-k+1)
      end do
    else
      zi(:) = xi(:)
      zg(:) = g(:)
    end if

    ! Merge coincident or non increasing coordinates
    call merge_coincident(zi,zg,npi)
    if ( npi == 1 ) then
      f(:) = zg(1)
      return
    end if

    do j = 1, npo
      if ( xo(j) < zi(1) ) then
        f(j) = zg(1) + ex1*(zg(1)-zg(2))/(zi(1)-zi(2)) * (xo(j)-zi(1))
        cycle
      else if ( xo(j) >= zi(npi) ) then
        f(j) = zg(npi) + ex2*(zg(npi)-zg(npi-1))/(zi(npi)-zi(npi-1)) * &
               (xo(j)-zi(npi))
        cycle
      end if
      if ( npi == 2 ) then
        f(j) = (zg(1)*(zi(2)-xo(j)) + zg(2)*(xo(j)-zi(1))) / &
               (zi(2)-zi(1))
        cycle
      end if
      ir = 0
      do jj = 1, npi
        if ( xo(j) >= zi(jj) ) ir = ir + 1
      end do
      if ( ir == 1 ) then
        fmm = 2.0_rk8 * zg(1) - zg(2)
        xmm = 2.0_rk8 * zi(1) - zi(2)
        fpp = zg(ir+2)
        xpp = zi(ir+2)
      else if ( ir == (npi-1) ) then
        fpp = 2.0_rk8 * zg(npi) - zg(npi-1)
        xpp = 2.0_rk8 * zi(npi) - zi(npi-1)
        fmm = zg(ir-1)
        xmm = zi(ir-1)
      else
        fmm = zg(ir-1)
        xmm = zi(ir-1)
        fpp = zg(ir+2)
        xpp = zi(ir+2)
      end if
      fm     = zg(ir)
      xm     = zi(ir)
      fp     = zg(ir+1)
      xp     = zi(ir+1)
      delx   = xp - xm
      delxp  = xpp - xp
      delxm  = xm - xmm
      delx1  = xo(j) - xm
      delx2  = xp - xo(j)
      delxs  = delx**2
      delx1s = delx1**2
      delx2s = delx2**2
      spl = fm*(delx2/delx + delx1*delx2s/(delxs*delxm) - delx1s*     &
            delx2/((delx+delxp)*delxs)) + fp*(delx1/delx +            &
            delx1s*delx2/(delxs*delxp) - delx1*delx2s/((delx+delxm)*  &
            delxs)) - fmm * delx1*delx2s/((delx+delxm)*delx*delxm) -  &
            fpp * delx1s*delx2/((delx+delxp)*delx*delxp)
      clin = (fm*delx2 + fp*delx1)/delx
      f(j) = alfa*clin + (1.0_rk8-alfa)*spl
    end do
  end subroutine interp1d_r8

  pure subroutine interp1d_r4_new(xi,g,xo,f,alfa,ex1,ex2,zi,zg,sxi,sxo)
    !$acc routine seq
    implicit none
    integer(ik4), intent(in) :: sxi, sxo
    real(rk4), dimension(sxi), intent(in) :: xi, g
    real(rk4), dimension(sxo), intent(in) :: xo
    real(rk4), dimension(sxo), intent(out) :: f
    real(rk4), dimension(sxi), intent(out) :: zi,zg
    real(rk4), intent(in) :: alfa
    real(rk4), intent(in) :: ex1
    real(rk4), intent(in) :: ex2

    real(rk4) :: fmm, fpp, xmm, xpp
    real(rk4) :: fm, xm, fp, xp, delx, delxp, delxm
    real(rk4) :: delx1, delx2, delxs, delx1s, delx2s
    real(rk4) :: spl, clin

    integer(ik4) :: npi, npo
    integer(ik4) :: k, j, jj, ir

    npi = sxi
    npo = sxo
    if ( xi(1) >= xi(npi) ) then
      do k = 1, npi
        zi(k) = xi(npi-k+1)
        zg(k) = g(npi-k+1)
      end do
    else
      zi(:) = xi(:)
      zg(:) = g(:)
    end if

    ! Merge coincident or non increasing coordinates
    call merge_coincident(zi,zg,npi)
    if ( npi == 1 ) then
      f(:) = zg(1)
      return
    end if

    do j = 1, npo
      if ( xo(j) < zi(1) ) then
        f(j) = zg(1) + ex1*(zg(1)-zg(2))/(zi(1)-zi(2)) * (xo(j)-zi(1))
        cycle
      else if ( xo(j) >= zi(npi) ) then
        f(j) = zg(npi) + ex2*(zg(npi)-zg(npi-1))/(zi(npi)-zi(npi-1)) * &
               (xo(j)-zi(npi))
        cycle
      end if
      if ( npi == 2 ) then
        f(j) = (zg(1)*(zi(2)-xo(j)) + zg(2)*(xo(j)-zi(1))) / &
               (zi(2)-zi(1))
        cycle
      end if
      ir = 0
      do jj = 1, npi
        if ( xo(j) >= zi(jj) ) ir = ir + 1
      end do
      if ( ir == 1 ) then
        fmm = 2.0_rk4 * zg(1) - zg(2)
        xmm = 2.0_rk4 * zi(1) - zi(2)
        fpp = zg(ir+2)
        xpp = zi(ir+2)
      else if ( ir == (npi-1) ) then
        fpp = 2.0_rk4 * zg(npi) - zg(npi-1)
        xpp = 2.0_rk4 * zi(npi) - zi(npi-1)
        fmm = zg(ir-1)
        xmm = zi(ir-1)
      else
        fmm = zg(ir-1)
        xmm = zi(ir-1)
        fpp = zg(ir+2)
        xpp = zi(ir+2)
      end if
      fm     = zg(ir)
      xm     = zi(ir)
      fp     = zg(ir+1)
      xp     = zi(ir+1)
      delx   = xp - xm
      delxp  = xpp - xp
      delxm  = xm - xmm
      delx1  = xo(j) - xm
      delx2  = xp - xo(j)
      delxs  = delx**2
      delx1s = delx1**2
      delx2s = delx2**2
      spl = fm*(delx2/delx + delx1*delx2s/(delxs*delxm) - delx1s*     &
            delx2/((delx+delxp)*delxs)) + fp*(delx1/delx +            &
            delx1s*delx2/(delxs*delxp) - delx1*delx2s/((delx+delxm)*  &
            delxs)) - fmm * delx1*delx2s/((delx+delxm)*delx*delxm) -  &
            fpp * delx1s*delx2/((delx+delxp)*delx*delxp)
      clin = (fm*delx2 + fp*delx1)/delx
      f(j) = alfa*clin + (1.0_rk4-alfa)*spl
    end do
  end subroutine interp1d_r4_new

  pure subroutine interp1d_r8_new(xi,g,xo,f,alfa,ex1,ex2,zi,zg,sxi,sxo)
    !$acc routine seq
    implicit none
    integer(ik4), intent(in) :: sxi, sxo
    real(rk8), dimension(sxi), intent(in) :: xi, g
    real(rk8), dimension(sxo), intent(in) :: xo
    real(rk8), dimension(sxo), intent(out) :: f
    real(rk8), dimension(sxi), intent(out) :: zi,zg
    real(rk8), intent(in) :: alfa
    real(rk8), intent(in) :: ex1
    real(rk8), intent(in) :: ex2

    real(rk8) :: fmm, fpp, xmm, xpp
    real(rk8) :: fm, xm, fp, xp, delx, delxp, delxm
    real(rk8) :: delx1, delx2, delxs, delx1s, delx2s
    real(rk8) :: spl, clin

    integer(ik4) :: npi, npo
    integer(ik4) :: k, j, jj, ir

    npi = sxi
    npo = sxo
    if ( xi(1) >= xi(npi) ) then
      do k = 1, npi
        zi(k) = xi(npi-k+1)
        zg(k) = g(npi-k+1)
      end do
    else
      zi(:) = xi(:)
      zg(:) = g(:)
    end if

    ! Merge coincident or non increasing coordinates
    call merge_coincident(zi,zg,npi)
    if ( npi == 1 ) then
      f(:) = zg(1)
      return
    end if

    do j = 1, npo
      if ( xo(j) < zi(1) ) then
        f(j) = zg(1) + ex1*(zg(1)-zg(2))/(zi(1)-zi(2)) * (xo(j)-zi(1))
        cycle
      else if ( xo(j) >= zi(npi) ) then
        f(j) = zg(npi) + ex2*(zg(npi)-zg(npi-1))/(zi(npi)-zi(npi-1)) * &
               (xo(j)-zi(npi))
        cycle
      end if
      if ( npi == 2 ) then
        f(j) = (zg(1)*(zi(2)-xo(j)) + zg(2)*(xo(j)-zi(1))) / &
               (zi(2)-zi(1))
        cycle
      end if
      ir = 0
      do jj = 1, npi
        if ( xo(j) >= zi(jj) ) ir = ir + 1
      end do
      if ( ir == 1 ) then
        fmm = 2.0_rk8 * zg(1) - zg(2)
        xmm = 2.0_rk8 * zi(1) - zi(2)
        fpp = zg(ir+2)
        xpp = zi(ir+2)
      else if ( ir == (npi-1) ) then
        fpp = 2.0_rk8 * zg(npi) - zg(npi-1)
        xpp = 2.0_rk8 * zi(npi) - zi(npi-1)
        fmm = zg(ir-1)
        xmm = zi(ir-1)
      else
        fmm = zg(ir-1)
        xmm = zi(ir-1)
        fpp = zg(ir+2)
        xpp = zi(ir+2)
      end if
      fm     = zg(ir)
      xm     = zi(ir)
      fp     = zg(ir+1)
      xp     = zi(ir+1)
      delx   = xp - xm
      delxp  = xpp - xp
      delxm  = xm - xmm
      delx1  = xo(j) - xm
      delx2  = xp - xo(j)
      delxs  = delx**2
      delx1s = delx1**2
      delx2s = delx2**2
      spl = fm*(delx2/delx + delx1*delx2s/(delxs*delxm) - delx1s*     &
            delx2/((delx+delxp)*delxs)) + fp*(delx1/delx +            &
            delx1s*delx2/(delxs*delxp) - delx1*delx2s/((delx+delxm)*  &
            delxs)) - fmm * delx1*delx2s/((delx+delxm)*delx*delxm) -  &
            fpp * delx1s*delx2/((delx+delxp)*delx*delxp)
      clin = (fm*delx2 + fp*delx1)/delx
      f(j) = alfa*clin + (1.0_rk8-alfa)*spl
    end do
  end subroutine interp1d_r8_new

  !
  ! Merge the points whose coordinate is not larger than the previous one
  ! into one point with the mean coordinate and value, until the
  ! coordinates are strictly increasing. npi is updated.
  !
  pure subroutine merge_coincident_k4(zi,zg,npi)
    !$acc routine seq
    implicit none
    real(rk4), dimension(:), intent(inout) :: zi, zg
    integer(ik4), intent(inout) :: npi
    integer(ik4) :: k, kk
    k = 2
    do while ( k <= npi )
      if ( zi(k) <= zi(k-1) ) then
        zi(k-1) = 0.5_rk4 * (zi(k-1) + zi(k))
        zg(k-1) = 0.5_rk4 * (zg(k-1) + zg(k))
        do kk = k, npi-1
          zi(kk) = zi(kk+1)
          zg(kk) = zg(kk+1)
        end do
        npi = npi - 1
        ! The merged point may now be out of order with the previous one
        k = max(k-1,2)
      else
        k = k + 1
      end if
    end do
  end subroutine merge_coincident_k4

  !
  ! Merge the points whose coordinate is not larger than the previous one
  ! into one point with the mean coordinate and value, until the
  ! coordinates are strictly increasing. npi is updated.
  !
  pure subroutine merge_coincident_k8(zi,zg,npi)
    !$acc routine seq
    implicit none
    real(rk8), dimension(:), intent(inout) :: zi, zg
    integer(ik4), intent(inout) :: npi
    integer(ik4) :: k, kk
    k = 2
    do while ( k <= npi )
      if ( zi(k) <= zi(k-1) ) then
        zi(k-1) = 0.5_rk8 * (zi(k-1) + zi(k))
        zg(k-1) = 0.5_rk8 * (zg(k-1) + zg(k))
        do kk = k, npi-1
          zi(kk) = zi(kk+1)
          zg(kk) = zg(kk+1)
        end do
        npi = npi - 1
        ! The merged point may now be out of order with the previous one
        k = max(k-1,2)
      else
        k = k + 1
      end if
    end do
  end subroutine merge_coincident_k8

end module mod_interp

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
