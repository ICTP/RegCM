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

module mod_zita

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_dynparam , only : mo_a0 , mo_h , mo_ztop

  implicit none

  private

  interface zita_interp
    module procedure zh3d
    module procedure zh4d
  end interface zita_interp

  public :: model_zitaf , model_zitah , zitasigma , sigmazita
  public :: md_zeta , md_zeta_h , md_fmz , md_fmz_h , gzita
  public :: md_ak , md_bk
  public :: zita_interp

  contains

  pure real(rkx) elemental function sigmazita(zita)
    implicit none
    real(rkx) , intent(in) :: zita
    sigmazita = 1.0_rkx - zita/mo_ztop
  end function sigmazita

  pure real(rkx) elemental function zitasigma(sigma)
    implicit none
    real(rkx) , intent(in) :: sigma
    zitasigma = mo_ztop * (1.0_rkx - sigma)
  end function zitasigma

  subroutine model_zitaf(zitaf)
    implicit none
    real(rkx) , intent(out) , dimension(:) :: zitaf
    real(rkx) :: dz
    integer :: kzp1 , kz , k
    kzp1 = size(zitaf)
    kz = kzp1-1
    dz = mo_ztop/real(kz,rkx)
    zitaf(kzp1) = 0.0_rkx
    do k = kz , 1 , -1
      zitaf(k) = zitaf(k+1) + dz
    end do
  end subroutine model_zitaf

  subroutine model_zitah(zitah)
    real(rkx) , intent(out) , dimension(:) :: zitah
    real(rkx) :: dz
    integer :: kz , k
    kz = size(zitah)
    dz = mo_ztop/real(kz,rkx)
    zitah(kz) = dz*0.5_rkx
    do k = kz-1 , 1 , -1
      zitah(k) = zitah(k+1) + dz
    end do
  end subroutine model_zitah

  ! Decay function

  pure real(rkx) elemental function zfz( )
    implicit none
    zfz = mo_ztop/(exp(mo_ztop/mo_h)-1.0_rkx)
  end function zfz

  pure real(rkx) elemental function bzita(zita)
    implicit none
    real(rkx) , intent(in) :: zita
    bzita = zfz( )*(exp(zita/mo_h)-1.0_rkx)
  end function bzita

  pure real(rkx) elemental function bzitap(zita)
    implicit none
    real(rkx) , intent(in) :: zita
    bzitap = zfz( )*exp(zita/mo_h)/mo_h
  end function bzitap

  pure real(rkx) elemental function gzita(zita)
    implicit none
    real(rkx) , intent(in) :: zita
    real(rkx) :: ratio
    ratio = zita/mo_ztop
    gzita = ((0.0_rkx - 1.0_rkx * mo_a0) * ratio**1 - &
             (3.0_rkx - 2.0_rkx * mo_a0) * ratio**2 + &
             (2.0_rkx - 1.0_rkx * mo_a0) * ratio**3) + 1.0_rkx
    !gzita = sinh((mo_ztop-zita)/mo_h)/sinh(mo_ztop/mo_h)
  end function gzita

  ! Derivative of decay function
  pure real(rkx) elemental function gzitap(zita)
    implicit none
    real(rkx) , intent(in) :: zita
    real(rkx) :: ratio
    ratio = zita/mo_ztop
    gzitap = ((0.0_rkx - 1.0_rkx * mo_a0) * ratio**0 - &
              (6.0_rkx - 4.0_rkx * mo_a0) * ratio**1 + &
              (6.0_rkx - 3.0_rkx * mo_a0) * ratio**2) / mo_ztop
    !gzitap = -(1.0_rkx/mo_h)*cosh((mo_ztop-zita)/mo_h)/sinh(mo_ztop/mo_h)
  end function gzitap

  ! Factor used to transform the vertical derivatives in zeta
  pure real(rkx) function md_fmz_h(zita,orog)
    implicit none
    real(rkx) , intent(in) :: zita , orog
    ! Equation 9
    md_fmz_h = 1.0_rkx/(gzitap(zita)*orog + bzitap(zita))
  end function md_fmz_h

  ! Elevation above orography as function of zita
  pure real(rkx) function md_zeta_h(zita,orog)
    implicit none
    real(rkx) , intent(in) :: zita , orog
    md_zeta_h = orog*(gzita(zita)-1.0_rkx) + bzita(zita)
  end function md_zeta_h

  pure real(rkx) elemental function md_ak(zita)
    implicit none
    real(rkx) , intent(in) :: zita
    md_ak = bzita(zita)
  end function md_ak

  pure real(rkx) elemental function md_bk(zita)
    implicit none
    real(rkx) , intent(in) :: zita
    md_bk = gzita(zita)
  end function md_bk

  pure real(rkx) function md_fmz(zita,geopot)
    implicit none
    real(rkx) , intent(in) :: zita , geopot
    md_fmz = md_fmz_h(zita,geopot*regrav)
  end function md_fmz

  ! Elevation above orography as function of zita
  pure real(rkx) function md_zeta(zita,geopot)
    implicit none
    real(rkx) , intent(in) :: zita , geopot
    md_zeta = md_zeta_h(zita,geopot*regrav)
  end function md_zeta

  subroutine zh3d(nx1,nx2,ny1,ny2,nz,f,zeta,tvirt,sigmah,ps,imet)
    implicit none
    integer(ik4) , intent(in) :: nx1 , nx2 , ny1 , ny2 , nz , imet
    real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: f
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: zeta
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: tvirt
    real(rkx) , dimension(:) , pointer , intent(in) :: sigmah
    real(rkx) , dimension(:,:) , pointer , intent(in) :: ps
    real(rkx) , dimension(nz) :: psigma , pz
    real(rkx) , dimension(nz) :: fz
    integer(ik4) :: i , j , k , kk , ik
    real(rkx) :: iw1 , iw2
    if ( imet == 1 ) then
      do j = ny1 , ny2
        do i = nx1 , nx2
          psigma = ps(i,j) * sigmah
          pz = ps(i,j)*exp(-egrav*zeta(i,j,:)/rgas/tvirt(i,j,:))
          do k = 1 , nz
            if ( pz(k) < psigma(1) ) then
              fz(k) = f(i,j,1)
            else if ( pz(k) > psigma(nz) ) then
              fz(k) = f(i,j,nz)
            else
              ! Find requested pressure level
              ik = 1
              do kk = 2 , nz
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
      do j = ny1 , ny2
        do i = nx1 , nx2
          psigma = ps(i,j) * sigmah
          pz = ps(i,j)*exp(-egrav*zeta(i,j,:)/rgas/tvirt(i,j,:))
          do k = 1 , nz
            if ( pz(k) < psigma(1) ) then
              fz(k) = f(i,j,1)
            else if ( pz(k) > psigma(nz) ) then
              fz(k) = f(i,j,nz)
            else
              ! Find requested pressure level
              ik = 1
              do kk = 2 , nz
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
    implicit none
    integer(ik4) , intent(in) :: nx1 , nx2 , ny1 , ny2 , nz , nn , imet
    real(rkx) , dimension(:,:,:,:) , pointer , intent(inout) :: f
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: zeta
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: tvirt
    real(rkx) , dimension(:) , pointer , intent(in) :: sigmah
    real(rkx) , dimension(:,:) , pointer , intent(in) :: ps
    real(rkx) , dimension(nz) :: psigma , pz
    real(rkx) , dimension(nz) :: fz
    integer(ik4) :: i , j , k , n , kk , ik
    real(rkx) :: iw1 , iw2
    if ( imet == 1 ) then
      do n = 1 , nn
        do j = ny1 , ny2
          do i = nx1 , nx2
            psigma = ps(i,j) * sigmah
            pz = ps(i,j)*exp(-egrav*zeta(i,j,:)/rgas/tvirt(i,j,:))
            do k = 1 , nz
              if ( pz(k) < psigma(1) ) then
                fz(k) = f(i,j,1,n)
              else if ( pz(k) > psigma(nz) ) then
                fz(k) = f(i,j,nz,n)
              else
                ! Find requested pressure level
                ik = 1
                do kk = 2 , nz
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
      do n = 1 , nn
        do j = ny1 , ny2
          do i = nx1 , nx2
            psigma = ps(i,j) * sigmah
            pz = ps(i,j)*exp(-egrav*zeta(i,j,:)/rgas/tvirt(i,j,:))
            do k = 1 , nz
              if ( pz(k) < psigma(1) ) then
                fz(k) = f(i,j,1,n)
              else if ( pz(k) > psigma(nz) ) then
                fz(k) = f(i,j,nz,n)
              else
                ! Find requested pressure level
                ik = 1
                do kk = 2 , nz
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
