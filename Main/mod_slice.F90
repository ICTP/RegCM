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
  use mod_memutil
  use mod_atm_interface
  use mod_che_interface
  use mod_pbl_interface
  use mod_rad_interface

  implicit none

  private

  public :: mkslice , init_slice

  integer(ik4) :: ix1 , ix2 , jx1 , jx2
  integer(ik4) :: id1 , id2 , jd1 , jd2
  real(rkx) , dimension(:,:) , pointer :: rpsb
  real(rkx) , dimension(:,:) , pointer :: rpsdotb
  contains

  subroutine init_slice
    implicit none
    if ( idiffu == 1 ) then
      ix1 = ice1gb
      ix2 = ice2gb
      jx1 = jce1gb
      jx2 = jce2gb
      id1 = ide1gb
      id2 = ide2gb
      jd1 = jde1gb
      jd2 = jde2gb
    else if ( idiffu == 2 ) then
      ix1 = ice1ga
      ix2 = ice2ga
      jx1 = jce1ga
      jx2 = jce2ga
      id1 = ide1ga
      id2 = ide2ga
      jd1 = jde1ga
      jd2 = jde2ga
    else if ( idiffu == 3 ) then
      ix1 = ice1gc
      ix2 = ice2gc
      jx1 = jce1gc
      jx2 = jce2gc
      id1 = ide1gc
      id2 = ide2gc
      jd1 = jde1gc
      jd2 = jde2gc
    end if
    call getmem2d(rpsb,jx1,jx2,ix1,ix2,'slice:rpsb')
    call getmem2d(rpsdotb,jd1,jd2,id1,id2,'slice:rpsdotb')
  end subroutine init_slice

  subroutine mkslice
    implicit none
    real(rkx) :: cell
    integer(ik4) :: i , j , k , n

    do concurrent ( j = jx1:jx2 , i = ix1:ix2 )
      rpsb(j,i) = d_one/sfs%psb(j,i)
    end do
    do concurrent ( j = jd1:jd2 , i = id1:id2 )
      rpsdotb(j,i) = d_one/sfs%psdotb(j,i)
    end do

    do concurrent ( j = jd1:jd2 , i = id1:id2 , k = 1:kz )
      atms%ubd3d(j,i,k) = atm2%u(j,i,k)*rpsdotb(j,i)
      atms%vbd3d(j,i,k) = atm2%v(j,i,k)*rpsdotb(j,i)
    end do

    do concurrent ( j = jce1:jce2 , i = ice1:ice2 , k = 1:kz )
      atms%ubx3d(j,i,k) = d_rfour *                        &
              (atms%ubd3d(j,i,k)   + atms%ubd3d(j,i+1,k) + &
               atms%ubd3d(j+1,i,k) + atms%ubd3d(j+1,i+1,k))
      atms%vbx3d(j,i,k) = d_rfour *                        &
              (atms%vbd3d(j,i,k)   + atms%vbd3d(j,i+1,k) + &
               atms%vbd3d(j+1,i,k) + atms%vbd3d(j+1,i+1,k))
    end do
    do concurrent ( j = jx1:jx2 , i = ix1:ix2 , k = 1:kz )
      atms%tb3d(j,i,k) = atm2%t(j,i,k)*rpsb(j,i)
      atms%qxb3d(j,i,k,iqv) = max(atm2%qx(j,i,k,iqv)*rpsb(j,i),minqq)
    end do
    do n = iqfrst , iqlst
      do k = 1 , kz
        do i = ix1 , ix2
          do j = jx1 , jx2
            atms%qxb3d(j,i,k,n) = atm2%qx(j,i,k,n)*rpsb(j,i)
            if ( atms%qxb3d(j,i,k,n) < minqq * minqq ) then
              atms%qxb3d(j,i,k,n) = d_zero
            end if
          end do
        end do
      end do
    end do
    if ( ichem == 1 ) then
      do n = 1 , ntr
        do k = 1 , kz
          do i = ix1 , ix2
            do j = jx1 , jx2
              atms%chib3d(j,i,k,n) = atm2%chi(j,i,k,n)*rpsb(j,i)
              if ( atms%chib3d(j,i,k,n) < mintr ) then
                atms%chib3d(j,i,k,n) = d_zero
              end if
            end do
          end do
        end do
      end do
    end if

    do concurrent ( j = jce1ga:jce2ga , i = ice1ga:ice2ga , k = 1:kz )
      atms%tv3d(j,i,k) = atms%tb3d(j,i,k) * &
              (d_one + ep1*atms%qxb3d(j,i,k,iqv) - atms%qxb3d(j,i,k,iqc))
    end do

    if ( idynamic == 2 ) then
      do concurrent ( j = jx1:jx2 , i = ix1:ix2 , k = 1:kz )
        atms%ppb3d(j,i,k) = atm2%pp(j,i,k)*rpsb(j,i)
      end do
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 , k = 2:kz )
        atms%pb3d(j,i,k) = atm0%pr(j,i,k) + atms%ppb3d(j,i,k)
      end do
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 )
        atms%pb3d(j,i,1) = max(atm0%pr(j,i,1) + atms%ppb3d(j,i,1), &
                          ptop*d_1000+1.0_rkx)
      end do
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 )
        atms%ps2d(j,i) = atm0%ps(j,i) + ptop*d_1000 + atms%ppb3d(j,i,kz)
      end do
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 )
        atms%pf3d(j,i,1) = ptop*d_1000
      end do
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 )
        atms%pf3d(j,i,kzp1) = atms%ps2d(j,i)
      end do
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 , k = 2:kz )
        atms%pf3d(j,i,k) = atm0%pf(j,i,k) + &
                   d_half*(atms%ppb3d(j,i,k-1)+atms%ppb3d(j,i,k))
      end do
    else
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 , k = 1:kz )
        atms%pb3d(j,i,k) = (hsigma(k)*sfs%psb(j,i) + ptop)*d_1000
      end do
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 )
        atms%ps2d(j,i) = (sfs%psb(j,i)+ptop)*d_1000
      end do
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 , k = 1:kzp1 )
        atms%pf3d(j,i,k) = (sigma(k)*sfs%psb(j,i) + ptop)*d_1000
      end do
    end if

    do concurrent ( j = jci1:jci2 , i = ici1:ici2 )
      atms%rhox2d(j,i) = atms%ps2d(j,i)/(rgas*atms%tb3d(j,i,kz))
    end do

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atms%rhob3d(j,i,k) = atms%pb3d(j,i,k)/(rgas*atms%tb3d(j,i,k))
          atms%th3d(j,i,k) = atms%tb3d(j,i,k) * &
                          (p00/atms%pb3d(j,i,k))**rovcp
          atms%tp3d(j,i,k) = atms%tb3d(j,i,k) * &
                          (atms%ps2d(j,i)/atms%pb3d(j,i,k))**rovcp
        end do
      end do
    end do

    if ( idynamic == 2 ) then
      do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz )
        atms%wpx3d(j,i,k) = omega(j,i,k)
      end do
      do concurrent ( j = jx1:jx2 , i = ix1:ix2 , k = 1:kzp1 )
        atms%wb3d(j,i,k) = atm2%w(j,i,k)*rpsb(j,i)
      end do
    else
      ! Omega in the hydrostatic model is in cb/s
      do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz )
        atms%wpx3d(j,i,k) = omega(j,i,k) * d_1000 ! Pa/s
      end do
      atms%wb3d(:,:,1) = d_zero
      atms%wb3d(:,:,kzp1) = d_zero
      do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 2:kz )
        atms%wb3d(j,i,k) = -d_half*regrav * &
                   (atms%wpx3d(j,i,k-1)/atms%rhob3d(j,i,k-1) + &
                    atms%wpx3d(j,i,k)/atms%rhob3d(j,i,k))
      end do
    end if

    !
    ! Find 700 mb theta
    !
    if ( icldmstrat == 1 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          atms%th700(j,i) = atms%th3d(j,i,kz)
          vertloop: &
          do k = 1 , kz-1
            if ( atms%pb3d(j,i,k) > 70000.0 ) then
              atms%th700(j,i) = twt(k,1) * atms%th3d(j,i,k+1) + &
                                twt(k,2) * atms%th3d(j,i,k)
              exit vertloop
            end if
          end do vertloop
        end do
      end do
    end if
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
    if ( idynamic == 1 ) then
      !
      ! compute the height at full (za) and half (zq) sigma levels:
      ! CONSTANT FOR NON-HYDROSTATIC
      !
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 )
        atms%zq(j,i,kzp1) = d_zero
      end do
      do concurrent ( j = jce1ga:jce2ga , i = ice1ga:ice2ga , k = kz:1:-1 )
        cell = ptop * rpsb(j,i)
        atms%zq(j,i,k) = atms%zq(j,i,k+1) + rovg * atms%tv3d(j,i,k) *  &
                      log((sigma(k+1)+cell)/(sigma(k)+cell))
      end do
      do concurrent ( j = jce1ga:jce2ga , i = ice1ga:ice2ga , k = 1:kz )
        atms%za(j,i,k) = d_half*(atms%zq(j,i,k) + atms%zq(j,i,k+1))
      end do
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 , k = 1:kz )
        atms%dzq(j,i,k) = atms%zq(j,i,k) - atms%zq(j,i,k+1)
      end do
    end if

    do concurrent ( j = jce1:jce2 , i = ice1:ice2 , k = 1:kz )
      atms%qsb3d(j,i,k) = pfwsat(atms%tb3d(j,i,k),atms%pb3d(j,i,k))
    end do
    do concurrent ( j = jce1:jce2 , i = ice1:ice2 , k = 1:kz )
      atms%rhb3d(j,i,k) = atms%qxb3d(j,i,k,iqv)/atms%qsb3d(j,i,k)
      atms%rhb3d(j,i,k) = min(max(atms%rhb3d(j,i,k),rhmin),rhmax)
    end do

    contains

#include <pfesat.inc>
#include <pfwsat.inc>

  end subroutine mkslice

end module mod_slice

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
