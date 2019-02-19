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
  real(rkx) , pointer , dimension(:,:,:) :: ux , vx , wx
  real(rkx) , pointer , dimension(:,:,:) :: tkex

  public allocate_moloch

  contains

#include <cpmf.inc>

  subroutine allocate_moloch
    implicit none
    call getmem3d(s,jci1,jci2,ici1,ici2,1,kzp1,'moloch:s')
    call getmem3d(ux,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'moloch:ux')
    call getmem3d(vx,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'moloch:vx')
    call getmem3d(wx,jci1,jci2,ici1,ici2,1,kz,'moloch:wx')
    if ( ibltyp == 2 ) then
      call getmem3d(tkex,jci1,jci2,ici1,ici2,1,kz,'moloch:tkex')
    end if
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

      subroutine advection
        implicit none
        integer(ik4) :: i , j , k , n
        real(rkx) :: dz
        real(rkx) , pointer , dimension(:,:,:) :: ptr

        dz = hzita / real(kz,rkx)

        call exchange_lr(mo_atm%u,2,jde1,jde2,ide1,ide2,1,kz)
        call exchange_bt(mo_atm%v,2,jde1,jde2,ide1,ide2,1,kz)

        ! Compute U-wind on T points

        do k = 1 , kz
          do i = ice1 , ice2
            do j = jcii1 , jci2
              ux(j,i,k) = 0.5625_rkx * (mo_atm%u(j,i,k)+mo_atm%u(j-1,i,k)) - &
                          0.0625_rkx * (mo_atm%u(j+1,i,k)+mo_atm%u(j-2,i,k))
            end do
          end do
        end do
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ice1 , ice2
              ux(jce1,i,k) = mo_atm%u(jde1,i,k)
              ux(jci1,i,k) = 0.5_rkx * (mo_atm%u(jde1,i,k)+mo_atm%u(jdi1,i,k))
            end do
          end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ice1 , ice2
              ux(jce2,i,k) = 0.5_rkx*(mo_atm%u(jde2,i,k) + mo_atm%u(jdi2,i,k))
            end do
          end do
        end if

        ! Compute V-wind on T points

        do k = 1 , kz
          do i = icii1 , jci2
            do j = jce1 , jce2
              vx(j,i,k) = 0.5625_rkx * (mo_atm%v(j,i,k)+mo_atm%v(j,i-1,k)) - &
                          0.0625_rkx * (mo_atm%v(j,i+1,k)+mo_atm%v(j,i-2,k))
            end do
          end do
        end do
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              vx(j,ice1,k) = mo_atm%v(j,ide1,k)
              vx(j,ici1,k) = 0.5_rkx * (mo_atm%v(j,ide1,k)+mo_atm%v(j,idi1,k))
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do i = ice1 , ice2
              vx(j,ice2,k) = 0.5_rkx*(mo_atm%v(j,ide2,k) + mo_atm%v(j,idi2,k))
            end do
          end do
        end if

        ! Compute W (and TKE if required) on zita levels

        do k = 2 , kz-1
          do i = ici1 , ici2
            do j = jci1 , jci2
              wx(j,i,k) = 0.5625_rkx * (mo_atm%w(j,i,k+1)+mo_atm%w(j,i,k)) - &
                          0.0625_rkx * (mo_atm%w(j,i,k+2)+mo_atm%w(j,i,k-1))
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            wx(j,i,1) = 0.5_rkx * (mo_atm%w(j,i,2)+mo_atm%w(j,i,1))
            wx(j,i,kz) = 0.5_rkx * (mo_atm%w(j,i,kzp1)+mo_atm%w(j,i,kz))
          end do
        end do

        if ( ibltyp == 2 ) then
          do k = 2 , kz-1
            do i = ici1 , ici2
              do j = jci1 , jci2
                tkex(j,i,k) = 0.5625_rkx * &
                              (mo_atm%tke(j,i,k+1)+mo_atm%tke(j,i,k)) - &
                              0.0625_rkx * &
                              (mo_atm%tke(j,i,k+2)+mo_atm%tke(j,i,k-1))
              end do
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              tkex(j,i,1) = 0.5_rkx*(mo_atm%tke(j,i,2)+mo_atm%tke(j,i,1))
              tkex(j,i,kz) = 0.5_rkx*(mo_atm%tke(j,i,kzp1)+mo_atm%tke(j,i,kz))
            end do
          end do
        end if

        call wafone(mo_atm%tetav,mo_atm%u,mo_atm%v,dx,dx,dz,dt, &
                    mddom%clv,mddom%fmyu)
        call wafone(mo_atm%pai,mo_atm%u,mo_atm%v,dx,dx,dz,dt, &
                    mddom%clv,mddom%fmyu)
        call wafone(ux,mo_atm%u,mo_atm%v,dx,dx,dz,dt, &
                    mddom%clv,mddom%fmyu)
        call wafone(vx,mo_atm%u,mo_atm%v,dx,dx,dz,dt, &
                    mddom%clv,mddom%fmyu)
        call wafone(wx,mo_atm%u,mo_atm%v,dx,dx,dz,dt, &
                    mddom%clv,mddom%fmyu)
        do n = 1 , nqx
          call assignpnt(mo_atm%qx,ptr,n)
          call wafone(ptr,mo_atm%u,mo_atm%v,dx,dx,dz,dt, &
                      mddom%clv,mddom%fmyu)
        end do
        if ( ibltyp == 2 ) then
          call wafone(tkex,mo_atm%u,mo_atm%v,dx,dx,dz,dt, &
                      mddom%clv,mddom%fmyu)
        end if
        if ( ichem == 1 ) then
          do n = 1 , ntr
            call assignpnt(mo_atm%trac,ptr,n)
            call wafone(ptr,mo_atm%u,mo_atm%v,dx,dx,dz,dt, &
                        mddom%clv,mddom%fmyu)
          end do
        end if

        call exchange_lr(ux,2,jce1,jce2,ice1,ice2,1,kz)
        call exchange_bt(vx,2,jce1,jce2,ice1,ice2,1,kz)

        ! Back to wind points: U (fourth order)

        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdii2
              mo_atm%u(j,i,k) = 0.5625_rkx * (ux(j,i,k)+ux(j+1,i,k)) - &
                                0.0625_rkx * (ux(j-1,i,k)+ux(j+2,i,k))
            end do
          end do
        end do
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = idi1 , idi2
              mo_atm%u(jdi2,i,k) = 0.5_rkx * (ux(jci2,i,k)+ux(jce2,i,k))
            end do
          end do
        end if

        ! Back to wind points: V (fourth order)

        do k = 1 , kz
          do i = idi1 , idii2
            do j = jdi1 , jdi2
              mo_atm%v(j,i,k) = 0.5625_rkx * (vx(j,i,k)+vx(j,i+1,k)) - &
                                0.0625_rkx * (vx(j,i-1,k)+vx(j,i+2,k))
            end do
          end do
        end do
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jdi1 , jdi2
              mo_atm%v(j,idi2,k) = 0.5_rkx * (vx(j,ici2,k)+vx(j,ice2,k))
            end do
          end do
        end if

        ! Back to half-levels

        do k = 3 , kz - 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              mo_atm%w(j,i,k) = 0.5625_rkx * (wx(j,i,k)+wx(j,i,k-1)) - &
                                0.0625_rkx * (wx(j,i,k+1)+wx(j,i,k-2))
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            mo_atm%w(j,i,2) = 0.5_rkx * (wx(j,i,2)+wx(j,i,1))
            mo_atm%w(j,i,kz) = 0.5_rkx * (wx(j,i,kz)+wx(j,i,kzm1))
          end do
        end do

        if ( ibltyp == 2 ) then
          do k = 3 , kz - 1
            do i = ici1 , ici2
              do j = jci1 , jci2
                mo_atm%tke(j,i,k) = 0.5625_rkx*(tkex(j,i,k)+tkex(j,i,k-1)) - &
                                    0.0625_rkx*(tkex(j,i,k+1)+tkex(j,i,k-2))
              end do
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              mo_atm%tke(j,i,2) = 0.5_rkx * (tkex(j,i,2)+tkex(j,i,1))
              mo_atm%tke(j,i,kz) = 0.5_rkx * (tkex(j,i,kz)+tkex(j,i,kzm1))
            end do
          end do
        end if
      end subroutine advection

      subroutine wafone(p,u,v,dx,dy,dz,dt,clv,fmyu)
        implicit none
        real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: p
        real(rkx) , dimension(:,:,:) , pointer , intent(in) :: u , v
        real(rkx) , dimension(:,:) , pointer , intent(in) :: clv , fmyu
        real(rkx) , intent(in) :: dx , dy , dz , dt
        integer(ik4) :: j , i , k
        integer(ik4) :: k1 , k1m1 , ih , ihm1 , im1 , jh , jhm1 , jm1
        real(rkx) , dimension(jci1:jci2,1:kzp1) :: wfw
        real(rkx) , dimension(jci1:jci2,ici1gb:ici2gb,1:kz) :: wz
        real(rkx) , dimension(jci1gb:jci2gb,ici1:ici2,1:kz) :: p0
        real(rkx) , dimension(jci1:jci2,ici1:ice2) :: zpby
        real(rkx) , dimension(jci1:jce2) :: zpbw
        real(rkx) :: zamu , zcost , r , b , zphi , is , zdv
        real(rkx) :: zcostx , zcosty , zhxvt , zhxvtn

        ! Vertical advection

        zcost = dt / dz
        do j = jci1 , jci2
          wfw(j,1) = d_zero
          wfw(j,kzp1) = d_zero
        end do

        do i = ici1 , ici2
          do k = 2 , kz
            do j = jci1 , jci2
              zamu = s(j,i,k)*zcost
              if ( zamu >= d_zero ) then
                is = d_one
                k1 = k - 1
                k1m1 = k1 - 1
                k1m1 = max(k1m1,1)
              else
                is = -d_one
                k1 = k + 1
                k1m1 = k1 - 1
                k1 = min(k1,kz)
              end if
              r = rdeno(p(j,i,k1),p(j,i,k1m1),p(j,i,k),p(j,i,k-1))
              b = max(d_zero, min(d_two, max(r, min(d_two*r,d_one))))
              zphi = is + zamu * b - is * b
              wfw(j,k) = d_half * zamu*((d_one+zphi)*p(j,i,k-1) + &
                                        (d_one-zphi)*p(j,i,k))
            end do
          end do
          do k = 1 , kz
            do j = jci1 , jci2
              zdv = (s(j,i,k+1)-s(j,i,k)) * zcost
              wz(j,i,k) = p(j,i,k) + wfw(j,k) - wfw(j,k+1) + p(j,i,k)*zdv
            end do
          end do
        end do

        call exchange_tb(wz,2,jci1,jci2,ici1,ici2,1,kz)

        ! Meridional advection

        zcosty = dt/dy

        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jce2
              zamu = v(j,i,k)*zcosty
              if ( zamu > d_zero ) then
                is = d_one
                ih = max(i-1,icross1+1)
              else
                is = -d_one
                ih = min(i+1,icross2-1)
              end if
              ihm1 = max(ih-1,icross1+1)
              im1 = max(i-1,icross1+1)
              r = rdeno(wz(j,ih,k), wz(j,ihm1,k), wz(j,i,k), wz(j,im1,k))
              b = max(d_zero, min(d_two, max(r, min(d_two*r,d_one))))
              zphi = is+zamu*b -is*b
              zpby(j,i) = d_half*zamu * &
                ((d_one+zphi)*wz(j,im1,k)+(d_one-zphi)*wz(j,i,k))
            end do
          end do
          do i = ici1 , ici2
            zhxvtn = clv(j,i+1)*fmyu(j,i)
            zhxvt  = clv(j,i)*fmyu(j,i)
            do j = jci1 , jci2
              zdv = (v(j,i+1,k)*zhxvtn - v(j,i,k)*zhxvt)*zcosty
              p0(j,i,k) = wz(j,i,k) + &
                      zpby(j,i)*zhxvt - zpby(j,i+1)*zhxvtn + p(j,i,k)*zdv
            end do
          end do
        end do

        call exchange_lr(p0,2,jci1,jci2,ici1,ici2,1,kz)

        ! Zonal advection

        do k = 1 , kz
          do i = ici1 , ici2
            zcostx = dt*fmyu(j,i)/dx
            do j = jci1 , jce2
              zamu = u(j-1,i,k)*zcostx
              if ( zamu > d_zero ) then
                is = d_one
                jh = max(j-1,jcross1+1)
              else
                is = -d_one
                jh = min(j+1,jcross2-1)
              end if
              jhm1 = max(jh-1,jcross1+1)
              jm1 = max(j-1,jcross1+1)
              r = rdeno (p0(jh,i,k), p0(jhm1,i,k), p0(j,i,k), p0(jm1,i,k))
              b = max(d_zero, min(d_two, max(r, min(d_two*r,d_one))))
              zphi = is+zamu*b -is*b
              zpbw(j) = d_half*zamu * &
                   ((d_one+zphi)*p0(j-1,i,k)+(d_one-zphi)*p0(j,i,k))
            end do
            do j = jci1 , jci2
              p(j,i,k) = p0(j,i,k) + zpbw(j) - zpbw(j+1) + &
                        p(j,i,k)*zdv
            end do
          end do
        end do
      end subroutine wafone

  end subroutine moloch

end module mod_moloch

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
