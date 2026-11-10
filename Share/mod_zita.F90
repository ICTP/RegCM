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

module mod_zita

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_dynparam, only : mo_lcap, mo_b0

  implicit none

  private

  interface zita_interp
    module procedure zh3d_r8
    module procedure zh3d_r4
    module procedure zh4d_r8
    module procedure zh4d_r4
  end interface zita_interp

  ! With two levels (CLM stand alone forcing only) the lowest layer is
  ! 80 m deep over flat terrain, so the lowest half level is near 40 m.
  real(rkx), parameter :: clmsa_dz = 80.0_rkx
  ! Uncapped coordinate: finite stand in for the infinite height of the
  ! top interface (zita = h), so that it can be written in the output.
  ! No computation is meant to use it.
  real(rkx), parameter :: uncapped_ztop = 100000.0_rkx

  public :: model_zitaf, model_zitah, model_dzitah, zitasigma, sigmazita
  public :: md_zeta, md_zeta_h, md_fmz, md_fmz_h, bzita, gzita
  public :: md_ak, md_bk
  public :: zita_interp

  contains

  pure real(rkx) elemental function sigmazita(zita,ztop)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: zita, ztop
    sigmazita = 1.0_rkx - zita/ztop
  end function sigmazita

  pure real(rkx) elemental function zitasigma(sigma,ztop)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: sigma, ztop
    zitasigma = ztop * (1.0_rkx - sigma)
  end function zitasigma

  subroutine model_zitaf(zitaf,ztop,zh)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: ztop, zh
    real(rkx), intent(out), dimension(:) :: zitaf
    real(rkx) :: dz
    integer :: kzp1, kz, k
    kzp1 = size(zitaf)
    kz = kzp1-1
    zitaf(kzp1) = 0.0_rkx
    zitaf(1) = ztop
    if ( kz == 2 ) then
      zitaf(2) = clmsa_zitaf(ztop,zh)
      return
    end if
    dz = ztop/real(kz,rkx)
    do k = kz, 2, -1
      zitaf(k) = zitaf(k+1) + dz
    end do
  end subroutine model_zitaf

  subroutine model_zitah(zitah,ztop,zh)
!$acc routine seq
    real(rkx), intent(in) :: ztop, zh
    real(rkx), intent(out), dimension(:) :: zitah
    real(rkx) :: dz
    integer :: kz, k
    kz = size(zitah)
    if ( kz == 2 ) then
      zitah(2) = clmsa_zitaf(ztop,zh)*0.5_rkx
      zitah(1) = (ztop + clmsa_zitaf(ztop,zh))*0.5_rkx
      return
    end if
    dz = ztop/real(kz,rkx)
    zitah(kz) = dz*0.5_rkx
    zitah(1) = ztop - dz*0.5_rkx
    do k = kz-1, 2, -1
      zitah(k) = zitah(k+1) + dz
    end do
  end subroutine model_zitah

  ! Zita of the interface at clmsa_dz above flat terrain (bzita inverse).
  ! The uncapped height function has no closed inverse: Newton iterations
  ! from the near surface value zita = clmsa_dz/b0.
  pure real(rkx) function clmsa_zitaf(ztop,zh)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: ztop, zh
    integer(ik4) :: it
    if ( mo_lcap ) then
      clmsa_zitaf = zh*log(1.0_rkx + clmsa_dz/zfz(ztop,zh))
    else
      clmsa_zitaf = clmsa_dz/mo_b0
      do it = 1, 20
        clmsa_zitaf = clmsa_zitaf - &
          (bzita(clmsa_zitaf,ztop,zh)-clmsa_dz)/bzitap(clmsa_zitaf,ztop,zh)
      end do
    end if
  end function clmsa_zitaf

  ! Distance in zita between the half levels k and k+1, for k < kz.
  pure real(rkx) function model_dzitah(k,kz,ztop)
!$acc routine seq
    implicit none
    integer(ik4), intent(in) :: k, kz
    real(rkx), intent(in) :: ztop
    if ( kz == 2 ) then
      model_dzitah = ztop*0.5_rkx
    else
      model_dzitah = ztop/real(kz,rkx)
    end if
    if ( k < 1 .or. k >= kz ) model_dzitah = 0.0_rkx
  end function model_dzitah

  ! Decay function

  pure real(rkx) elemental function zfz(ztop,zh)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: ztop, zh
    zfz = ztop/(exp(ztop/zh)-1.0_rkx)
  end function zfz

  ! Height above flat terrain as function of zita. Capped: exponential
  ! stretching from 0 to ztop. Uncapped (original MOLOCH, ztop = h):
  ! z = -h*b(zita)*log(1-zita/h), b(zita) = b0+(1-b0)*zita/h, which is a
  ! reference pressure coordinate p = p0*(1-zita/h) for b0 = 1 and goes to
  ! infinity at zita = h.
  pure real(rkx) elemental function bzita(zita,ztop,zh)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: zita, ztop, zh
    real(rkx) :: ratio
    if ( mo_lcap ) then
      bzita = zfz(ztop,zh)*(exp(zita/zh)-1.0_rkx)
    else
      ratio = zita/ztop
      if ( ratio < 1.0_rkx ) then
        bzita = -ztop*(mo_b0+(1.0_rkx-mo_b0)*ratio)*log(1.0_rkx-ratio)
      else
        bzita = uncapped_ztop
      end if
    end if
  end function bzita

  ! Derivative of bzita in zita
  pure real(rkx) elemental function bzitap(zita,ztop,zh)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: zita, ztop, zh
    real(rkx) :: ratio
    if ( mo_lcap ) then
      bzitap = zfz(ztop,zh)*exp(zita/zh)/zh
    else
      ratio = zita/ztop
      if ( ratio < 1.0_rkx ) then
        bzitap = (mo_b0+(1.0_rkx-mo_b0)*ratio)/(1.0_rkx-ratio) - &
                 (1.0_rkx-mo_b0)*log(1.0_rkx-ratio)
      else
        ! Infinite at the top interface: md_fmz_h returns zero there
        bzitap = 0.0_rkx
      end if
    end if
  end function bzitap

  ! The decay function is the same in both coordinates (in the uncapped
  ! one ztop = h, as in the original MOLOCH).
  pure real(rkx) elemental function gzita(zita,ztop,a0)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: zita, ztop, a0
    real(rkx) :: ratio
    ratio = zita/ztop
    gzita = ((0.0_rkx - 1.0_rkx * a0) * ratio**1 - &
             (3.0_rkx - 2.0_rkx * a0) * ratio**2 + &
             (2.0_rkx - 1.0_rkx * a0) * ratio**3) + 1.0_rkx
    !gzita = sinh((ztop-zita)/zh)/sinh(ztop/zh)
  end function gzita

  ! Derivative of decay function
  pure real(rkx) elemental function gzitap(zita,ztop,a0)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: zita, ztop, a0
    real(rkx) :: ratio
    ratio = zita/ztop
    gzitap = ((0.0_rkx - 1.0_rkx * a0) * ratio**0 - &
              (6.0_rkx - 4.0_rkx * a0) * ratio**1 + &
              (6.0_rkx - 3.0_rkx * a0) * ratio**2) / ztop
    !gzitap = -(1.0_rkx/zh)*cosh((ztop-zita)/zh)/sinh(ztop/zh)
  end function gzitap

  ! Factor used to transform the vertical derivatives in zeta
  pure real(rkx) function md_fmz_h(zita,orog,ztop,zh,a0)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: zita, orog, ztop, zh, a0
    ! Equation 9
    if ( .not. mo_lcap .and. zita >= ztop ) then
      ! Uncapped top interface: infinite dz/dzita (as fmzh in MOLOCH)
      md_fmz_h = 0.0_rkx
    else
      md_fmz_h = 1.0_rkx/(gzitap(zita,ztop,a0)*orog + bzitap(zita,ztop,zh))
    end if
  end function md_fmz_h

  ! Elevation above orography as function of zita
  pure real(rkx) function md_zeta_h(zita,orog,ztop,zh,a0)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: zita, orog, ztop, zh, a0
    md_zeta_h = orog*(gzita(zita,ztop,a0)-1.0_rkx) + bzita(zita,ztop,zh)
  end function md_zeta_h

  pure real(rkx) elemental function md_ak(zita,ztop,zh)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: zita, ztop, zh
    md_ak = bzita(zita,ztop,zh)
  end function md_ak

  pure real(rkx) elemental function md_bk(zita,ztop,a0)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: zita, ztop, a0
    md_bk = gzita(zita,ztop,a0)
  end function md_bk

  pure real(rkx) function md_fmz(zita,geopot,ztop,zh,a0)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: zita, geopot, ztop, zh, a0
    md_fmz = md_fmz_h(zita,geopot*regrav,ztop,zh,a0)
  end function md_fmz

  ! Elevation above orography as function of zita
  pure real(rkx) function md_zeta(zita,geopot,ztop,zh,a0)
!$acc routine seq
    implicit none
    real(rkx), intent(in) :: zita, geopot, ztop, zh, a0
    md_zeta = md_zeta_h(zita,geopot*regrav,ztop,zh,a0)
  end function md_zeta

  subroutine zh3d_r8(nx1,nx2,ny1,ny2,nz,f,zeta,tvirt,sigmah,ps,imet)
    implicit none
    integer(ik4), intent(in) :: nx1, nx2, ny1, ny2, nz, imet
    real(rk8), dimension(:,:,:), pointer, contiguous, intent(inout) :: f
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: zeta
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: tvirt
    real(rkx), dimension(:), pointer, contiguous, intent(in) :: sigmah
    real(rkx), dimension(:,:), pointer, contiguous, intent(in) :: ps
    real(rk8), dimension(nz) :: psigma, pz
    real(rk8), dimension(nz) :: fz
    integer(ik4) :: i, j, k, kk, ik
    real(rk8) :: iw1, iw2
    if ( imet == 1 ) then
      do j = ny1, ny2
        do i = nx1, nx2
          psigma = ps(i,j) * sigmah
          pz = ps(i,j)*exp(-egrav*zeta(i,j,:)/rgas/tvirt(i,j,:))
          do k = 1, nz
            if ( pz(k) < psigma(1) ) then
              fz(k) = f(i,j,1)
            else if ( pz(k) > psigma(nz) ) then
              fz(k) = f(i,j,nz)
            else
              ! Find requested pressure level
              ik = 1
              do kk = 2, nz
                if ( pz(k) < psigma(kk) ) then
                  ik = kk
                  exit
                end if
              end do
              iw1 = log(pz(k)/psigma(ik-1))/log(psigma(ik)/psigma(ik-1))
              iw2 = 1.0_rkx - iw1
              fz(k) = iw1 * f(i,j,ik) + iw2 * f(i,j,ik-1)
            end if
          end do
          do k = 1, nz
            f(i,j,k) = fz(k)
          end do
        end do
      end do
    else
      do j = ny1, ny2
        do i = nx1, nx2
          psigma = ps(i,j) * sigmah
          pz = ps(i,j)*exp(-egrav*zeta(i,j,:)/rgas/tvirt(i,j,:))
          do k = 1, nz
            if ( pz(k) < psigma(1) ) then
              fz(k) = f(i,j,1)
            else if ( pz(k) > psigma(nz) ) then
              fz(k) = f(i,j,nz)
            else
              ! Find requested pressure level
              ik = 1
              do kk = 2, nz
                if ( pz(k) < psigma(kk) ) then
                  ik = kk
                  exit
                end if
              end do
              iw1 = (pz(k)-psigma(ik-1))/(psigma(ik)-psigma(ik-1))
              iw2 = 1.0_rkx - iw1
              fz(k) = iw1 * f(i,j,ik) + iw2 * f(i,j,ik-1)
            end if
          end do
          do k = 1, nz
            f(i,j,k) = fz(k)
          end do
        end do
      end do
    end if
  end subroutine zh3d_r8

  subroutine zh3d_r4(nx1,nx2,ny1,ny2,nz,f,zeta,tvirt,sigmah,ps,imet)
    implicit none
    integer(ik4), intent(in) :: nx1, nx2, ny1, ny2, nz, imet
    real(rk4), dimension(:,:,:), pointer, contiguous, intent(inout) :: f
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: zeta
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: tvirt
    real(rkx), dimension(:), pointer, contiguous, intent(in) :: sigmah
    real(rkx), dimension(:,:), pointer, contiguous, intent(in) :: ps
    real(rk8), dimension(nz) :: psigma, pz
    real(rk8), dimension(nz) :: fz
    integer(ik4) :: i, j, k, kk, ik
    real(rk8) :: iw1, iw2
    if ( imet == 1 ) then
      do j = ny1, ny2
        do i = nx1, nx2
          psigma = ps(i,j) * sigmah
          pz = ps(i,j)*exp(-egrav*zeta(i,j,:)/rgas/tvirt(i,j,:))
          do k = 1, nz
            if ( pz(k) < psigma(1) ) then
              fz(k) = f(i,j,1)
            else if ( pz(k) > psigma(nz) ) then
              fz(k) = f(i,j,nz)
            else
              ! Find requested pressure level
              ik = 1
              do kk = 2, nz
                if ( pz(k) < psigma(kk) ) then
                  ik = kk
                  exit
                end if
              end do
              iw1 = log(pz(k)/psigma(ik-1))/log(psigma(ik)/psigma(ik-1))
              iw2 = 1.0_rkx - iw1
              fz(k) = iw1 * f(i,j,ik) + iw2 * f(i,j,ik-1)
            end if
          end do
          do k = 1, nz
            f(i,j,k) = real(fz(k),rk4)
          end do
        end do
      end do
    else
      do j = ny1, ny2
        do i = nx1, nx2
          psigma = ps(i,j) * sigmah
          pz = ps(i,j)*exp(-egrav*zeta(i,j,:)/rgas/tvirt(i,j,:))
          do k = 1, nz
            if ( pz(k) < psigma(1) ) then
              fz(k) = f(i,j,1)
            else if ( pz(k) > psigma(nz) ) then
              fz(k) = f(i,j,nz)
            else
              ! Find requested pressure level
              ik = 1
              do kk = 2, nz
                if ( pz(k) < psigma(kk) ) then
                  ik = kk
                  exit
                end if
              end do
              iw1 = (pz(k)-psigma(ik-1))/(psigma(ik)-psigma(ik-1))
              iw2 = 1.0_rkx - iw1
              fz(k) = iw1 * f(i,j,ik) + iw2 * f(i,j,ik-1)
            end if
          end do
          do k = 1, nz
            f(i,j,k) = real(fz(k),rk4)
          end do
        end do
      end do
    end if
  end subroutine zh3d_r4

  subroutine zh4d_r8(nx1,nx2,ny1,ny2,nz,nn,f,zeta,tvirt,sigmah,ps,imet)
    implicit none
    integer(ik4), intent(in) :: nx1, nx2, ny1, ny2, nz, nn, imet
    real(rk8), dimension(:,:,:,:), pointer, contiguous, intent(inout) :: f
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: zeta
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: tvirt
    real(rkx), dimension(:), pointer, contiguous, intent(in) :: sigmah
    real(rkx), dimension(:,:), pointer, contiguous, intent(in) :: ps
    real(rk8), dimension(nz) :: psigma, pz
    real(rk8), dimension(nz) :: fz
    integer(ik4) :: i, j, k, n, kk, ik
    real(rk8) :: iw1, iw2
    if ( imet == 1 ) then
      do n = 1, nn
        do j = ny1, ny2
          do i = nx1, nx2
            psigma = ps(i,j) * sigmah
            pz = ps(i,j)*exp(-egrav*zeta(i,j,:)/rgas/tvirt(i,j,:))
            do k = 1, nz
              if ( pz(k) < psigma(1) ) then
                fz(k) = f(i,j,1,n)
              else if ( pz(k) > psigma(nz) ) then
                fz(k) = f(i,j,nz,n)
              else
                ! Find requested pressure level
                ik = 1
                do kk = 2, nz
                  if ( pz(k) < psigma(kk) ) then
                    ik = kk
                    exit
                  end if
                end do
                iw1 = log(pz(k)/psigma(ik-1))/log(psigma(ik)/psigma(ik-1))
                iw2 = 1.0_rkx - iw1
                fz(k) = iw1 * f(i,j,ik,n) + iw2 * f(i,j,ik-1,n)
              end if
            end do
            do k = 1, nz
              f(i,j,k,n) = fz(k)
            end do
          end do
        end do
      end do
    else
      do n = 1, nn
        do j = ny1, ny2
          do i = nx1, nx2
            psigma = ps(i,j) * sigmah
            pz = ps(i,j)*exp(-egrav*zeta(i,j,:)/rgas/tvirt(i,j,:))
            do k = 1, nz
              if ( pz(k) < psigma(1) ) then
                fz(k) = f(i,j,1,n)
              else if ( pz(k) > psigma(nz) ) then
                fz(k) = f(i,j,nz,n)
              else
                ! Find requested pressure level
                ik = 1
                do kk = 2, nz
                  if ( pz(k) < psigma(kk) ) then
                    ik = kk
                    exit
                  end if
                end do
                iw1 = (pz(k)-psigma(ik-1))/(psigma(ik)-psigma(ik-1))
                iw2 = 1.0_rkx - iw1
                fz(k) = iw1 * f(i,j,ik,n) + iw2 * f(i,j,ik-1,n)
              end if
            end do
            do k = 1, nz
              f(i,j,k,n) = fz(k)
            end do
          end do
        end do
      end do
    end if
  end subroutine zh4d_r8

  subroutine zh4d_r4(nx1,nx2,ny1,ny2,nz,nn,f,zeta,tvirt,sigmah,ps,imet)
    implicit none
    integer(ik4), intent(in) :: nx1, nx2, ny1, ny2, nz, nn, imet
    real(rk4), dimension(:,:,:,:), pointer, contiguous, intent(inout) :: f
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: zeta
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: tvirt
    real(rkx), dimension(:), pointer, contiguous, intent(in) :: sigmah
    real(rkx), dimension(:,:), pointer, contiguous, intent(in) :: ps
    real(rk8), dimension(nz) :: psigma, pz
    real(rk8), dimension(nz) :: fz
    integer(ik4) :: i, j, k, n, kk, ik
    real(rk8) :: iw1, iw2
    if ( imet == 1 ) then
      do n = 1, nn
        do j = ny1, ny2
          do i = nx1, nx2
            psigma = ps(i,j) * sigmah
            pz = ps(i,j)*exp(-egrav*zeta(i,j,:)/rgas/tvirt(i,j,:))
            do k = 1, nz
              if ( pz(k) < psigma(1) ) then
                fz(k) = f(i,j,1,n)
              else if ( pz(k) > psigma(nz) ) then
                fz(k) = f(i,j,nz,n)
              else
                ! Find requested pressure level
                ik = 1
                do kk = 2, nz
                  if ( pz(k) < psigma(kk) ) then
                    ik = kk
                    exit
                  end if
                end do
                iw1 = log(pz(k)/psigma(ik-1))/log(psigma(ik)/psigma(ik-1))
                iw2 = 1.0_rkx - iw1
                fz(k) = iw1 * f(i,j,ik,n) + iw2 * f(i,j,ik-1,n)
              end if
            end do
            do k = 1, nz
              f(i,j,k,n) = real(fz(k),rk4)
            end do
          end do
        end do
      end do
    else
      do n = 1, nn
        do j = ny1, ny2
          do i = nx1, nx2
            psigma = ps(i,j) * sigmah
            pz = ps(i,j)*exp(-egrav*zeta(i,j,:)/rgas/tvirt(i,j,:))
            do k = 1, nz
              if ( pz(k) < psigma(1) ) then
                fz(k) = f(i,j,1,n)
              else if ( pz(k) > psigma(nz) ) then
                fz(k) = f(i,j,nz,n)
              else
                ! Find requested pressure level
                ik = 1
                do kk = 2, nz
                  if ( pz(k) < psigma(kk) ) then
                    ik = kk
                    exit
                  end if
                end do
                iw1 = (pz(k)-psigma(ik-1))/(psigma(ik)-psigma(ik-1))
                iw2 = 1.0_rkx - iw1
                fz(k) = iw1 * f(i,j,ik,n) + iw2 * f(i,j,ik-1,n)
              end if
            end do
            do k = 1, nz
              f(i,j,k,n) = real(fz(k),rk4)
            end do
          end do
        end do
      end do
    end if
  end subroutine zh4d_r4

end module mod_zita

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
