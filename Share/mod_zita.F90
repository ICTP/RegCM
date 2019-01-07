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

  implicit none

  private

  real(rkx) , parameter :: t0 = 280.0_rkx
  real(rkx) , parameter :: b0 = 0.5_rkx
  real(rkx) , parameter :: hzita = rgas*t0/egrav

  public :: rdeno , md_zeta , md_zeta_h , md_fmz , md_fmz_h , hzita
  public :: bzita , gzita
  public :: zita_interp

  contains

  ! WAF Scheme
  real(rkx) elemental function rdeno(t1, t2, t3, t4)
    implicit none
    real(rkx) , intent(in) :: t1 , t2 , t3 , t4
    real(rkx) :: zzden
    zzden = t3 - t4
    rdeno = (t1-t2) / sign(max(abs(zzden),1.e-15_rkx),zzden)
  end function rdeno

  ! Decay function
  real(rkx) elemental function gzita(zita)
    implicit none
    real(rkx) , intent(in) :: zita
    gzita = 1.0_rkx - 3.0_rkx*(zita/hzita)**2 + 2.0_rkx*(zita/hzita)**3
  end function gzita

  ! Derivative of decay function
  real(rkx) elemental function gzitap(zita)
    implicit none
    real(rkx) , intent(in) :: zita
    gzitap = (-6.0_rkx*(zita/hzita) + 6.0_rkx*(zita/hzita)**2)/hzita
  end function gzitap

  ! Stretching function
  real(rkx) elemental function bzita(zita)
    implicit none
    real(rkx) , intent(in) :: zita
    bzita = b0 + (1.0_rkx-b0)*(zita/hzita)
  end function bzita
  
  ! Derivative of stretching function
  real(rkx) elemental function bzitap(zita)
    implicit none
    real(rkx) , intent(in) :: zita
    bzitap = (1.0_rkx-b0)/hzita
  end function bzitap

  ! Factor used to transform the vertical derivatives in zeta
  real(rkx) function md_fmz(zita,geopot)
    implicit none
    real(rkx) :: zita , geopot
    real(rkx) :: zfz
    zfz = 1.0_rkx - zita/hzita
    ! Equation 9
    md_fmz = zfz /( bzita(zita) + geopot/egrav*zfz*gzitap(zita) - &
                    hzita*zfz*log(zfz)*bzitap(zita) )
  end function md_fmz

  real(rkx) function md_fmz_h(zita,orog)
    implicit none
    real(rkx) :: zita , orog
    real(rkx) :: zfz
    zfz = 1.0_rkx - zita/hzita
    ! Equation 9
    md_fmz_h = zfz /( bzita(zita) + orog*zfz*gzitap(zita) - &
                      hzita*zfz*log(zfz)*bzitap(zita) )
  end function md_fmz_h

  ! Elevation above orography as function of zita
  real(rkx) function md_zeta(zita,geopot)
    implicit none
    real(rkx) :: zita , geopot
    real(rkx) :: zfz , orog
    zfz = 1.0_rkx - zita/hzita
    orog = geopot/egrav
    ! Equation 7 with removal of orography and check for negatives
    md_zeta = max((orog*gzita(zita)-hzita*bzita(zita)*log(zfz))-orog, 0.0_rkx)
  end function md_zeta

  ! Elevation above orography as function of zita
  real(rkx) function md_zeta_h(zita,orog)
    implicit none
    real(rkx) :: zita , orog
    real(rkx) :: zfz
    zfz = 1.0_rkx - zita/hzita
    ! Equation 7 with removal of orography and check for negatives
    md_zeta_h = max((orog*gzita(zita)-hzita*bzita(zita)*log(zfz))-orog, 0.0_rkx)
  end function md_zeta_h

  subroutine zita_interp(nx,ny,nz,f,zeta,tvirt,sigmah,ps,ptop)
    implicit none
    integer(ik4) , intent(in) :: nx , ny , nz
    real(rkx) , dimension(:,:,:) , intent(inout) :: f
    real(rkx) , dimension(:,:,:) , intent(in) :: zeta
    real(rkx) , dimension(:,:,:) , intent(in) :: tvirt
    real(rkx) , dimension(:) , intent(in) :: sigmah
    real(rkx) , dimension(:,:) , intent(in) :: ps
    real(rkx) , intent(in) :: ptop
    real(rkx) , dimension(nz) :: psigma , pz
    real(rkx) , dimension(nz) :: fz
    integer(ik4) :: i , j , k , kk , ik
    real(rkx) :: iw1 , iw2
    do j = 1 , ny
      do i = 1 , nx
        psigma = (ps(i,j)-ptop)* sigmah + ptop
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
  end subroutine zita_interp

end module mod_zita

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
