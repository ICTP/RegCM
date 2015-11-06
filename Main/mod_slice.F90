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

module mod_slice
  !
  ! Fill 3D spaces for calculations
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_runparams
  use mod_atm_interface
  use mod_che_interface
  use mod_pbl_interface
  use mod_rad_interface

  implicit none

  private

  public :: mkslice

  contains

  subroutine mkslice
    implicit none
    real(rk8) :: cell
    integer(ik4) :: i , j , k , n
    real(rk8) , dimension(jce1:jce2,ice1:ice2) :: rpsb
    real(rk8) , dimension(jde1:jde2,ide1:ide2) :: rpsdotb

    do i = ice1 , ice2
      do j = jce1 , jce2
        rpsb(j,i) = d_one/sfs%psb(j,i)
      end do
    end do
    do i = ide1 , ide2
      do j = jde1 , jde2
        rpsdotb(j,i) = d_one/sfs%psdotb(j,i)
      end do
    end do

    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          atms%ubd3d(j,i,k) = atm2%u(j,i,k)*rpsdotb(j,i)
          atms%vbd3d(j,i,k) = atm2%v(j,i,k)*rpsdotb(j,i)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%ubx3d(j,i,k) = d_rfour*             &
              (atm2%u(j,i,k)   + atm2%u(j,i+1,k) + &
               atm2%u(j+1,i,k) + atm2%u(j+1,i+1,k)) * rpsb(j,i)
          atms%vbx3d(j,i,k) = d_rfour*             &
              (atm2%v(j,i,k)   + atm2%v(j,i+1,k) + &
               atm2%v(j+1,i,k) + atm2%v(j+1,i+1,k)) * rpsb(j,i)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%tb3d(j,i,k) = atm2%t(j,i,k)*rpsb(j,i)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%qxb3d(j,i,k,iqv) = max(atm2%qx(j,i,k,iqv),minqv)*rpsb(j,i)
        end do
      end do
    end do
    do n = iqc , nqx
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            if ( atm2%qx(j,i,k,n) > minqx ) then
              atms%qxb3d(j,i,k,n) = atm2%qx(j,i,k,n)*rpsb(j,i)
            else
              atms%qxb3d(j,i,k,n) = d_zero
            end if
          end do
        end do
      end do
    end do
    if ( ichem == 1 ) then
      do n = 1 , ntr
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              if ( chib(j,i,k,n) > mintr ) then
                atms%chib3d(j,i,k,n) = chib(j,i,k,n)*rpsb(j,i)
              else
                atms%chib3d(j,i,k,n) = d_zero
              end if
            end do
          end do
        end do
      end do
    end if

    if ( idynamic == 2 ) then
      do k = 1 , kzp1
        do i = ice1 , ice2
          do j = jce1 , jce2
            if ( abs(atm2%w(j,i,k)) > minww ) then
              atms%wb3d(j,i,k) = atm2%w(j,i,k)*rpsb(j,i)
            else
              atms%wb3d(j,i,k) = d_zero
            end if
          end do
        end do
      end do
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            atms%ppb3d(j,i,k) = atm2%pp(j,i,k)*rpsb(j,i)
          end do
        end do
      end do
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%ps2d(j,i) = atm0%ps(j,i) + ptop*d_1000 + atms%ppb3d(j,i,kz)
        end do
      end do
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%pf3d(j,i,1) = ptop*d_1000
          atms%pf3d(j,i,kzp1) = atms%ps2d(j,i)
        end do
      end do
      do k = 2 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            atms%pf3d(j,i,k) = atm0%pf(j,i,k) + &
                     d_half*(atms%ppb3d(j,i,k-1)+atms%ppb3d(j,i,k))
          end do
        end do
      end do
    else
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%ps2d(j,i) = (sfs%psb(j,i)+ptop)*d_1000
        end do
      end do
      do k = 1 , kzp1
        do i = ice1 , ice2
          do j = jce1 , jce2
            atms%pf3d(j,i,k) = (sigma(k)*sfs%psb(j,i) + ptop)*d_1000
          end do
        end do
      end do
    end if

    do i = ice1 , ice2
      do j = jce1 , jce2
        atms%rhox2d(j,i) = atms%ps2d(j,i)/(rgas*atms%tb3d(j,i,kz))
      end do
    end do

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%pb3d(j,i,k) = atm2%pr(j,i,k)
          atms%rhob3d(j,i,k) = atm2%rho(j,i,k)
          atms%th3d(j,i,k) = atms%tb3d(j,i,k) * &
                      (1.0D5/atms%pb3d(j,i,k))**rovcp
          atms%tp3d(j,i,k) = atms%tb3d(j,i,k) * &
                      (atms%ps2d(j,i)/atms%pb3d(j,i,k))**rovcp
        end do
      end do
    end do
    !
    ! Find tropopause hgt.
    !
    ktrop(:,:) = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        do k = kz , 1 , -1
          if ( atms%pb3d(j,i,k) < ptrop(j,i) ) then
            ktrop(j,i) = k
            exit
          end if
        end do
      end do
    end do
    !
    ! compute the height at full (za) and half (zq) sigma levels:
    !
    do i = ice1 , ice2
      do j = jce1 , jce2
        atms%zq(j,i,kzp1) = d_zero
      end do
    end do
    if ( idynamic == 2 ) then
      do k = kz , 1 , -1
        do i = ice1 , ice2
          do j = jce1 , jce2
            cell = ptop * rpsb(j,i)
            atms%zq(j,i,k) = atms%zq(j,i,k+1) + rovg * atm0%t(j,i,k) *  &
                      log((sigma(k+1)+cell)/(sigma(k)+cell))
          end do
        end do
      end do
      do i = ice1 , ice2
        do j = jce1 , jce2
          cell = ptop * rpsb(j,i)
          atms%za(j,i,kz) = rovg * atm0%t(j,i,kz) * &
                   log((sigma(kzp1)+cell)/(hsigma(kz)+cell))
        end do
      end do
      do k = kz-1 , 1 , -1
        do i = ice1 , ice2
          do j = jce1 , jce2
            cell = ptop * rpsb(j,i)
            atms%za(j,i,k) = atms%za(j,i,k+1) + rovg * atm0%t(j,i,k) * &
                     log((hsigma(k+1)+cell)/(hsigma(k)+cell))
          end do
        end do
      end do
    else
      do k = kz , 1 , -1
        do i = ice1 , ice2
          do j = jce1 , jce2
            cell = ptop * rpsb(j,i)
            atms%zq(j,i,k) = atms%zq(j,i,k+1) + rovg * atms%tb3d(j,i,k) *  &
                      log((sigma(k+1)+cell)/(sigma(k)+cell))
          end do
        end do
      end do
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            atms%za(j,i,k) = d_half*(atms%zq(j,i,k) + atms%zq(j,i,k+1))
          end do
        end do
      end do
    end if
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%dzq(j,i,k) = atms%zq(j,i,k) - atms%zq(j,i,k+1)
        end do
      end do
    end do

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%qsb3d(j,i,k) = pfqsat(atms%tb3d(j,i,k),atms%pb3d(j,i,k))
          atms%rhb3d(j,i,k) = atms%qxb3d(j,i,k,iqv)/atms%qsb3d(j,i,k)
          atms%rhb3d(j,i,k) = min(max(atms%rhb3d(j,i,k),1.D-3),rhmax)
        end do
      end do
    end do

    if ( ipptls == 2 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            atms%wpx3d(j,i,k) = omega(j,i,k)
          end do
        end do
      end do
    else
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            atms%wpx3d(j,i,k) = omega(j,i,k) * d_1000 ! Pa/s
          end do
        end do
      end do
    end if

    if ( ibltyp == 2 ) then
      do k = 1 , kzp1
        do i = ice1 , ice2
          do j = jce1 , jce2
            atms%tkeb3d(j,i,k) = atm2%tke(j,i,k)*rpsb(j,i)
          end do
        end do
      end do
    end if

  end subroutine mkslice

end module mod_slice

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
