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
  use mod_dynparam, only : mo_a0, mo_h, mo_ztop

  implicit none (type, external)

  private

  interface zita_interp
    module procedure zh3d
    module procedure zh4d
  end interface zita_interp

  public :: model_zitaf, model_zitah, zitasigma, sigmazita
  public :: md_zeta, md_zeta_h, md_fmz, md_fmz_h, gzita
  public :: md_ak, md_bk
  public :: zita_interp

  contains

  pure real(rkx) elemental function sigmazita(zita,ztop)
!$acc routine seq
    implicit none (type, external)
    real(rkx), intent(in) :: zita, ztop
    sigmazita = 1.0_rkx - zita/ztop
  end function sigmazita

  pure real(rkx) elemental function zitasigma(sigma,ztop)
!$acc routine seq
    implicit none (type, external)
    real(rkx), intent(in) :: sigma, ztop
    zitasigma = ztop * (1.0_rkx - sigma)
  end function zitasigma

  subroutine model_zitaf(zitaf,ztop)
!$acc routine seq
    implicit none (type, external)
    real(rkx), intent(in) :: ztop
    real(rkx), intent(out), dimension(:) :: zitaf
    real(rkx) :: dz
    integer :: kzp1, kz, k
    kzp1 = size(zitaf)
    kz = kzp1-1
    dz = ztop/real(kz,rkx)
    zitaf(kzp1) = 0.0_rkx
    do k = kz, 1, -1
      zitaf(k) = zitaf(k+1) + dz
    end do
  end subroutine model_zitaf

  subroutine model_zitah(zitah,ztop)
!$acc routine seq
    real(rkx), intent(in) :: ztop
    real(rkx), intent(out), dimension(:) :: zitah
    real(rkx) :: dz
    integer :: kz, k
    kz = size(zitah)
    dz = ztop/real(kz,rkx)
    zitah(kz) = dz*0.5_rkx
    do k = kz-1, 1, -1
      zitah(k) = zitah(k+1) + dz
    end do
  end subroutine model_zitah

  ! Decay function

  pure real(rkx) elemental function zfz(ztop,zh)
!$acc routine seq
    implicit none (type, external)
    real(rkx), intent(in) :: ztop, zh
    zfz = ztop/(exp(ztop/zh)-1.0_rkx)
  end function zfz

  pure real(rkx) elemental function bzita(zita,ztop,zh)
!$acc routine seq
    implicit none (type, external)
    real(rkx), intent(in) :: zita, ztop, zh
    bzita = zfz(ztop,zh)*(exp(zita/zh)-1.0_rkx)
  end function bzita

  pure real(rkx) elemental function bzitap(zita,ztop,zh)
!$acc routine seq
    implicit none (type, external)
    real(rkx), intent(in) :: zita, ztop, zh
    bzitap = zfz(ztop,zh)*exp(zita/zh)/zh
  end function bzitap

  pure real(rkx) elemental function gzita(zita,ztop,a0)
!$acc routine seq
    implicit none (type, external)
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
    implicit none (type, external)
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
    implicit none (type, external)
    real(rkx), intent(in) :: zita, orog, ztop, zh, a0
    ! Equation 9
    md_fmz_h = 1.0_rkx/(gzitap(zita,ztop,a0)*orog + bzitap(zita,ztop,zh))
  end function md_fmz_h

  ! Elevation above orography as function of zita
  pure real(rkx) function md_zeta_h(zita,orog,ztop,zh,a0)
!$acc routine seq
    implicit none (type, external)
    real(rkx), intent(in) :: zita, orog, ztop, zh, a0
    md_zeta_h = orog*(gzita(zita,ztop,a0)-1.0_rkx) + bzita(zita,ztop,zh)
  end function md_zeta_h

  pure real(rkx) elemental function md_ak(zita,ztop,zh)
!$acc routine seq
    implicit none (type, external)
    real(rkx), intent(in) :: zita, ztop, zh
    md_ak = bzita(zita,ztop,zh)
  end function md_ak

  pure real(rkx) elemental function md_bk(zita,ztop,a0)
!$acc routine seq
    implicit none (type, external)
    real(rkx), intent(in) :: zita, ztop, a0
    md_bk = gzita(zita,ztop,a0)
  end function md_bk

  pure real(rkx) function md_fmz(zita,geopot,ztop,zh,a0)
!$acc routine seq
    implicit none (type, external)
    real(rkx), intent(in) :: zita, geopot, ztop, zh, a0
    md_fmz = md_fmz_h(zita,geopot*regrav,ztop,zh,a0)
  end function md_fmz

  ! Elevation above orography as function of zita
  pure real(rkx) function md_zeta(zita,geopot,ztop,zh,a0)
!$acc routine seq
    implicit none (type, external)
    real(rkx), intent(in) :: zita, geopot, ztop, zh, a0
    md_zeta = md_zeta_h(zita,geopot*regrav,ztop,zh,a0)
  end function md_zeta

  subroutine zh3d(nx1,nx2,ny1,ny2,nz,f,zeta,tvirt,sigmah,ps,imet)
    implicit none (type, external)
    integer(ik4), intent(in) :: nx1, nx2, ny1, ny2, nz, imet
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(inout) :: f
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: zeta
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: tvirt
    real(rkx), dimension(:), pointer, contiguous, intent(in) :: sigmah
    real(rkx), dimension(:,:), pointer, contiguous, intent(in) :: ps
    real(rkx), dimension(nz) :: psigma, pz
    real(rkx), dimension(nz) :: fz
    integer(ik4) :: i, j, k, kk, ik
    real(rkx) :: iw1, iw2
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
          f(i,j,:) = fz
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
              iw1 = (pz(k)/psigma(ik-1))/(psigma(ik)/psigma(ik-1))
              iw2 = 1.0_rkx - iw1
              fz(k) = iw1 * f(i,j,ik) + iw2 * f(i,j,ik-1)
            end if
          end do
          f(i,j,:) = fz
        end do
      end do
    end if
  end subroutine zh3d

  subroutine zh4d(nx1,nx2,ny1,ny2,nz,nn,f,zeta,tvirt,sigmah,ps,imet)
    implicit none (type, external)
    integer(ik4), intent(in) :: nx1, nx2, ny1, ny2, nz, nn, imet
    real(rkx), dimension(:,:,:,:), pointer, contiguous, intent(inout) :: f
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: zeta
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(in) :: tvirt
    real(rkx), dimension(:), pointer, contiguous, intent(in) :: sigmah
    real(rkx), dimension(:,:), pointer, contiguous, intent(in) :: ps
    real(rkx), dimension(nz) :: psigma, pz
    real(rkx), dimension(nz) :: fz
    integer(ik4) :: i, j, k, n, kk, ik
    real(rkx) :: iw1, iw2
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
            f(i,j,:,n) = fz
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
                iw1 = (pz(k)/psigma(ik-1))/(psigma(ik)/psigma(ik-1))
                iw2 = 1.0_rkx - iw1
                fz(k) = iw1 * f(i,j,ik,n) + iw2 * f(i,j,ik-1,n)
              end if
            end do
            f(i,j,:,n) = fz
          end do
        end do
      end do
    end if
  end subroutine zh4d

end module mod_zita

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
