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
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_moloch

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_stdio
  use mod_service
  use mod_memutil
  use mod_atm_interface
  use mod_che_interface
  use mod_cu_interface
  use mod_lm_interface
  use mod_rad_interface
  use mod_pbl_interface
  use mod_micro_interface
  use mod_bdycod
  use mod_slice
  use mod_sun
  use mod_advection
  use mod_diffusion
  use mod_domain
  use mod_sladvection
  use mod_slabocean
  use mod_sound
  use mod_split
  use mod_timefilter
  use mod_massck
  use mod_zita

  implicit none

  private

  real(rkx) , pointer , dimension(:,:,:) :: s

  contains

#include <cpmf.inc>

  subroutine allocate_moloch
    implicit none
    call getmem3d(s,jci1,jci2,ici1,ici2,1,kzp1,'moloch:s')
  end subroutine allocate_moloch

  !
  ! Moloch dynamical integration engine
  !
  subroutine moloch
    implicit none
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'moloch'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

      subroutine filt2d(p,anu2,j1,j2,i1,i2)
        implicit none
        real(rkx) , dimension(:,:) , pointer , intent(inout) :: p
        real(rkx) , intent(in) :: anu2
        integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
        integer(ik4) :: j , i
        real(rkx) , dimension(j1:j2,i1:i2) :: p2
        call exchange(p,1,j1,j2,i1,i2)
        do i = i1 , i2
          do j = j1 , j2
            p2(j,i) = 0.125_rkx * (p(j,i-1)+p(j-1,i)+p(j+1,i)+p(j,i+1)) - &
                0.5_rkx*p(j,i)
          end do
        end do
        do i = i1 , i2
          do j = j1 , j2
            p(j,i) = p(j,i) + anu2 * p2(j,i)
          end do
        end do
      end subroutine filt2d

      subroutine filt3d(p,anu2,i1,i2,j1,j2)
        implicit none
        real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: p
        real(rkx) , intent(in) :: anu2
        integer(ik4) , intent(in) :: i1 , i2 , j1 , j2
        integer(ik4) :: j , i , k
        real(rkx) , dimension(j1:j2,i1:i2) :: p2
        call exchange(p,1,j1,j2,i1,i2,1,kz)
        do k = 1  , kz
          do i = i1 , i2
            do j = j1 , j2
              p2(j,i) = 0.125_rkx * (p(j,i-1,k) + p(j-1,i,k) + &
                                     p(j+1,i,k) + p(j,i+1,k)) - &
                        0.5_rkx   * p(j,i,k)
            end do
          end do
          do i = i1 , i2
            do j = j1 , j2
              p(j,i,k) = p(j,i,k) + anu2 * p2(j,i)
            end do
          end do
        end do
      end subroutine filt3d

      subroutine wafone(p,dz,dt,clv,j1,j2,i1,i2)
        implicit none
        real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: p
        real(rkx) , dimension(:,:) , pointer , intent(in) :: clv , fmyu
        real(rkx) , intent(in) :: dz , dt
        integer(ik4) , intent(in) :: i1 , i2 , j1 , j2
        integer(ik4) :: j , i , k
        integer(ik4) :: is , k1 , k1m1 , ih
        real(rkx) , dimension(j1:j2,1:kzp1) :: wfw
        real(rkx) , dimension(j1-1:j2+1,i1-1:i2+1,1:kz) :: wz
        real(rkx) , dimension(j1-1:j2+1,i1-1:i2+1,1:kz) :: p0
        real(rkx) , dimension(j1:j2,i1:i2) :: zpby
        real(rkx) :: zamu , zcost , r , b , zphi , is , zdv
        real(rkx) :: zcosty , zhxvtn , zhxvtn

        ! Vertical advection

        zcost = dt / dz
        do j = j1 , j2
          wfw(j,1) = d_zero
          wfw(j,kzp1) = d_zero
        end do

        do i = i1 , i2
          do k = 2 , kz
            do j = j1 , j2
              zamu = s(j,i,k)*zcost
              if ( zamu >= d_zero ) then
                is = d_one
                k1 = k - 1
                k1m1 = j1 - 1
                if ( k1m1 < 1 ) k1m1 = 1
              else
                is = -d_One
                k1 = k + 1
                k1m1 = k1 - 1
                if ( k1 > kz ) k1 = kz
              end if
              r = rdeno(p(j,i,k1),p(j,i,j1m1),p(j,i,k),p(j,i,k-1))
              b = max(d_zero, min(d_two, max(r, min(d_two*r,d_one))))
              zphi = is + zamu*b - is * b
              wfw(j,k) = d_half * zamu*((d_one+zphi)*p(j,i,k-1) + &
                                        (d_one-zphi)*p(j,i,k))
            end do
          end do
          do k = 1 , kz
            do j = j1 , j2
              zdv = (s(j,i,k+1)-s(j,i,k)) * zcost
              wz(j,i,k) = p(j,i,k) + wfw(j,k) - wfw(j,k+1) + p(j,i,k)*zdv
            end do
          end do
        end do
                
        call exchange(wz,1,j1,j2,i1,i2,1,kz)

        ! Meridional advection

        zcosty = dy/dy
        
        do k = 1 , kz
          do i = i1 , i2
            do j = j1 , j2
              zamu = mo_atm%v(j,i,k)*zcosty
              if ( zamu > d_zero ) then
                is = d_one
                ih = i - 1
              else
                is = - d_one
                ih = i + 1
              end if
              r = rdeno(wz(j,ih,k), wz(j,ih-1,k), wz(j,i,k), wz(j,i-1,k))
              b = max(d_zero, min(d_two, max(r, min(d_two*r,d_one)))) 
              zphi = is+zamu*b -is*b
              zpby(j,i) = d_half*zamu * &
                ((d_one+zphi)*wz(j,i-1,k)+(d_one-zphi)*wz(j,i,k))
            end do
          end do
          do i = i1 , i2
            zhxvtn = clv(j,i+1)*fmyu(j,i)
            zhxvt  = clv(j,i)*fmyu(j,i)
            do j = j1 , j2
              zdv = (mo_atm%v(j,i+1,k)*zhxvtn - mo_atm%v(j,i,k)*zhxvt)*zcosty
              p0(j,i,k) = wz(j,i,k) + &
                      zpby(j,i)*zhxvt -zpby(j,i+1)*zhxvtn + p(j,i,k)*zdv
            end do
          end do
        end do

        call exchange(p0,1,j1,j2,i1,i2,1,kz)

        ! Zonal advection

        do k = 1 , kz
          do i = i1 , i2
            zcostx = dt*fmyu(jlat)/dx
            do j = j1 , j2
              
            end do
          end do
        end do
      end subroutine wafone

  end subroutine moloch

end module mod_moloch

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
