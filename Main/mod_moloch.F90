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
  use mod_timefilter
  use mod_massck
  use mod_zita

  implicit none

  private

  real(rkx) , pointer , dimension(:,:,:) :: s
  real(rkx) , pointer , dimension(:,:,:) :: wx
  real(rkx) , pointer , dimension(:,:,:) :: tkex
  real(rkx) , pointer , dimension(:,:,:) :: wz
  real(rkx) , pointer , dimension(:,:,:) :: p0
  real(rkx) , pointer , dimension(:,:) :: zpby
  real(rkx) , pointer , dimension(:,:) :: zpbw

  real(rkx) , pointer , dimension(:,:,:) :: ten0
  real(rkx) , pointer , dimension(:,:,:) :: qen0
  real(rkx) , pointer , dimension(:,:,:,:) :: chiten0

  public :: allocate_moloch , moloch
  public :: uvstagtox , wstagtox

  contains

#include <cpmf.inc>

  subroutine allocate_moloch
    implicit none
    call getmem3d(s,jci1,jci2,ici1,ici2,1,kzp1,'moloch:s')
    call getmem3d(wx,jci1,jci2,ici1,ici2,1,kz,'moloch:wx')
    call getmem3d(wz,jci1,jci2,ici1gb,ici2gb,1,kz,'moloch:wz')
    call getmem3d(p0,jce1gb,jce2gb,ici1,ici2,1,kz,'moloch:p0')
    call getmem2d(zpby,jci1,jci2,ici1ga,ice2ga,'moloch:zpby')
    call getmem2d(zpbw,jci1ga,jce2ga,ici1,ici2,'moloch:zpbw')
    if ( ibltyp == 2 ) then
      call getmem3d(tkex,jci1,jci2,ici1,ici2,1,kz,'moloch:tkex')
    end if
    if ( idiag > 0 ) then
      call getmem3d(ten0,jci1,jci2,ici1,ici2,1,kz,'moloch:ten0')
      call getmem3d(qen0,jci1,jci2,ici1,ici2,1,kz,'moloch:qen0')
    end if
    if ( ichem == 1 ) then
      if ( ichdiag > 0 ) then
        call getmem4d(chiten0,jci1,jci2,ici1,ici2,1,kz,1,ntr,'moloch:chiten0')
      end if
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

    call boundary

    call advection

    !
    ! Next timestep ready : increment elapsed forecast time
    !
    call rcmtimer%advance( )
    if ( islab_ocean == 1 ) xslabtime = xslabtime + dtsec
    if ( rcmtimer%lcount == 2 ) then
      dtbat = dtsrf
      dt = dt2
      rdt = d_one/dt
      dtsq = dt*dt
      dtcb = dt*dt*dt
    end if
    !
    ! calculate new solar zenith angle
    !
    call zenitm(coszrs)

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

      subroutine boundary
        implicit none
        call exchange(mo_atm%u,1,jde1,jde2,ice1,ice2,1,kz)
        call exchange(mo_atm%v,1,jce1,jce2,ide1,ide2,1,kz)
        call exchange(mo_atm%t,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange(mo_atm%qx,1,jce1,jce2,ice1,ice2,1,kz,iqv,iqv)
        if ( ichem == 1 ) then
          call exchange(mo_atm%trac,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        end if

        if ( iboudy == 1 .or. iboudy == 5 ) then
          if ( idiag > 0 ) then
            ten0 = mo_atm%t(jci1:jci2,ici1:ici2,:)
            qen0 = mo_atm%qx(jci1:jci2,ici1:ici2,:,iqv)
          end if
          call nudge(iboudy,mo_atm%t,xtb)
          call nudge(iboudy,mo_atm%qx,xqb,iqv)
          call nudge(iboudy,mo_atm%u,mo_atm%v,xub,xvb)
          if ( idiag > 0 ) then
            tdiag%bdy = mo_atm%t(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%bdy = mo_atm%qx(jci1:jci2,ici1:ici2,:,iqv) - qen0
          end if
        else if ( iboudy == 4 ) then
          call sponge(mo_atm%t,xtb)
          call sponge(mo_atm%qx,xqb,iqv)
          call sponge(mo_atm%u,mo_atm%v,xub,xvb)
          if ( idiag > 0 ) then
            tdiag%bdy = mo_atm%t(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%bdy = mo_atm%qx(jci1:jci2,ici1:ici2,:,iqv) - qen0
          end if
        end if
        if ( ichem == 1 ) then
          if ( ichdiag > 0 ) then
            chiten0 = mo_atm%trac(jci1:jci2,ici1:ici2,:,:)
          end if
          if ( iboudy == 1 .or. iboudy == 5 ) then
            ! call nudge_chi(kz,mo_atm%trac)
          end if
          if ( ichdiag > 0 ) then
            cbdydiag = mo_atm%trac(jci1:jci2,ici1:ici2,:,:) - chiten0
          end if
        end if
      end subroutine boundary

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
        real(rkx) , pointer , dimension(:,:,:) :: ptr

        call uvstagtox(mo_atm%u,mo_atm%v,mo_atm%ux,mo_atm%vx)

        ! Compute W (and TKE if required) on zita levels

        call wstagtox(mo_atm%w,wx)

        if ( ibltyp == 2 ) then
          call wstagtox(mo_atm%tke,tkex)
        end if

        call wafone(mo_atm%tetav,mo_atm%u,mo_atm%v,dx,dx,mo_dz,dt, &
                    mddom%clv,mddom%fmyu)
        call wafone(mo_atm%pai,mo_atm%u,mo_atm%v,dx,dx,mo_dz,dt, &
                    mddom%clv,mddom%fmyu)
        call wafone(mo_atm%ux,mo_atm%u,mo_atm%v,dx,dx,mo_dz,dt, &
                    mddom%clv,mddom%fmyu)
        call wafone(mo_atm%vx,mo_atm%u,mo_atm%v,dx,dx,mo_dz,dt, &
                    mddom%clv,mddom%fmyu)
        call wafone(wx,mo_atm%u,mo_atm%v,dx,dx,mo_dz,dt, &
                    mddom%clv,mddom%fmyu)
        do n = 1 , nqx
          call assignpnt(mo_atm%qx,ptr,n)
          call wafone(ptr,mo_atm%u,mo_atm%v,dx,dx,mo_dz,dt, &
                      mddom%clv,mddom%fmyu)
        end do
        if ( ibltyp == 2 ) then
          call wafone(tkex,mo_atm%u,mo_atm%v,dx,dx,mo_dz,dt, &
                      mddom%clv,mddom%fmyu)
        end if
        if ( ichem == 1 ) then
          do n = 1 , ntr
            call assignpnt(mo_atm%trac,ptr,n)
            call wafone(ptr,mo_atm%u,mo_atm%v,dx,dx,mo_dz,dt, &
                        mddom%clv,mddom%fmyu)
          end do
        end if

        call exchange_lr(mo_atm%ux,2,jce1,jce2,ice1,ice2,1,kz)
        call exchange_bt(mo_atm%vx,2,jce1,jce2,ice1,ice2,1,kz)

        ! Back to wind points: U (fourth order)

        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jdi1 , jdii2-1
                mo_atm%u(j,i,k) = &
                  0.5625_rkx * (mo_atm%ux(j,i,k)+mo_atm%ux(j+1,i,k)) - &
                  0.0625_rkx * (mo_atm%ux(j-1,i,k)+mo_atm%ux(j+2,i,k))
              end do
            end do
          end do
          do k = 1 , kz
            do i = ici1 , ici2
              mo_atm%u(jdii2,i,k) = &
                0.5625_rkx * (mo_atm%ux(jci2,i,k)+mo_atm%ux(jce2,i,k)) - &
                0.0625_rkx * (mo_atm%ux(jcii2,i,k)+mo_atm%ux(jce2,i,k))
              mo_atm%u(jdi2,i,k) = &
                0.5_rkx * (mo_atm%ux(jci2,i,k)+mo_atm%ux(jce2,i,k))
            end do
          end do
        else
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jdi1 , jdii2
                mo_atm%u(j,i,k) = &
                  0.5625_rkx * (mo_atm%ux(j,i,k)+mo_atm%ux(j+1,i,k)) - &
                  0.0625_rkx * (mo_atm%ux(j-1,i,k)+mo_atm%ux(j+2,i,k))
              end do
            end do
          end do
        end if

        ! Back to wind points: V (fourth order)

        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do i = idi1 , idii2-1
              do j = jci1 , jci2
                mo_atm%v(j,i,k) = &
                  0.5625_rkx * (mo_atm%vx(j,i,k)+mo_atm%vx(j,i+1,k)) - &
                  0.0625_rkx * (mo_atm%vx(j,i-1,k)+mo_atm%vx(j,i+2,k))
              end do
            end do
          end do
          do k = 1 , kz
            do j = jci1 , jci2
              mo_atm%v(j,idii2,k) = &
                0.5625_rkx * (mo_atm%vx(j,ici2,k)+mo_atm%vx(j,ice2,k)) - &
                0.0625_rkx * (mo_atm%vx(j,icii2,k)+mo_atm%vx(j,ice2,k))
              mo_atm%v(j,idi2,k) = &
                0.5_rkx * (mo_atm%vx(j,ici2,k)+mo_atm%vx(j,ice2,k))
            end do
          end do
        else
          do k = 1 , kz
            do i = idi1 , idii2
              do j = jci1 , jci2
                mo_atm%v(j,i,k) = &
                  0.5625_rkx * (mo_atm%vx(j,i,k)+mo_atm%vx(j,i+1,k)) - &
                  0.0625_rkx * (mo_atm%vx(j,i-1,k)+mo_atm%vx(j,i+2,k))
              end do
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

        call exchange_bt(wz,2,jci1,jci2,ici1,ici2,1,kz)

        ! Meridional advection

        zcosty = dt/dy

        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
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

          call exchange_bt(zpby,1,jci1,jci2,ici1,ice2)

          do i = ici1 , ici2
            do j = jci1 , jci2
              zhxvtn = clv(j,i+1)*fmyu(j,i)
              zhxvt  = clv(j,i)*fmyu(j,i)
              zdv = (v(j,i+1,k)*zhxvtn - v(j,i,k)*zhxvt)*zcosty
              p0(j,i,k) = wz(j,i,k) + &
                      zpby(j,i)*zhxvt - zpby(j,i+1)*zhxvtn + p(j,i,k)*zdv
            end do
          end do
        end do

        call exchange_lr(p0,2,jci1,jci2,ici1,ici2,1,kz)

        if ( ma%has_bdyleft ) then
          p0(jce1,:,:) = p0(jci1,:,:)
        end if
        if ( ma%has_bdyright ) then
          p0(jce2,:,:) = p0(jci2,:,:)
        end if

        ! Zonal advection

        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jce2
              zcostx = dt*fmyu(j,i)/dx
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
              zpbw(j,i) = d_half*zamu * &
                   ((d_one+zphi)*p0(j-1,i,k)+(d_one-zphi)*p0(j,i,k))
            end do
          end do

          call exchange_lr(zpbw,1,jci1,jce2,ici1,ici2)

          do i = ici1 , ici2
            do j = jci1 , jci2
              p(j,i,k) = p0(j,i,k) + zpbw(j,i) - zpbw(j+1,i) + &
                        p(j,i,k)*zdv
            end do
          end do
        end do

      end subroutine wafone

  end subroutine moloch

  subroutine wstagtox(w,wx)
    implicit none
    real(rkx) , intent(in) , dimension(:,:,:) , pointer :: w
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: wx
    integer(ik4) :: i , j , k
    do k = 2 , kz-1
      do i = ici1 , ici2
        do j = jci1 , jci2
          wx(j,i,k) = 0.5625_rkx * (w(j,i,k+1)+w(j,i,k)) - &
                      0.0625_rkx * (w(j,i,k+2)+w(j,i,k-1))
        end do
      end do
    end do
    do i = ici1 , ici2
      do j = jci1 , jci2
        wx(j,i,1)  = 0.5_rkx * (w(j,i,2)+w(j,i,1))
        wx(j,i,kz) = 0.5_rkx * (w(j,i,kzp1)+w(j,i,kz))
      end do
    end do
  end subroutine wstagtox

  subroutine uvstagtox(u,v,ux,vx)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: u , v
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: ux , vx
    integer(ik4) :: i , j , k
    call exchange_lr(u,2,jde1,jde2,ice1,ice2,1,kz)
    call exchange_bt(v,2,jce1,jce2,ide1,ide2,1,kz)
    ! Compute U-wind on T points
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jcii1 , jci2
          ux(j,i,k) = 0.5625_rkx * (u(j,i,k)+u(j-1,i,k)) - &
                      0.0625_rkx * (u(j+1,i,k)+u(j-2,i,k))
        end do
      end do
    end do
    if ( ma%has_bdyleft ) then
      do k = 1 , kz
        do i = ice1 , ice2
          ux(jce1,i,k) = u(jde1,i,k)
          ux(jci1,i,k) = 0.5_rkx * (u(jde1,i,k)+u(jdi1,i,k))
        end do
      end do
    end if
    if ( ma%has_bdyright ) then
      do k = 1 , kz
        do i = ice1 , ice2
          ux(jce2,i,k) = 0.5_rkx*(u(jde2,i,k) + u(jdi2,i,k))
        end do
      end do
    end if
    ! Compute V-wind on T points
    do k = 1 , kz
      do i = icii1 , ici2
        do j = jce1 , jce2
          vx(j,i,k) = 0.5625_rkx * (v(j,i,k)+v(j,i-1,k)) - &
                      0.0625_rkx * (v(j,i+1,k)+v(j,i-2,k))
        end do
      end do
    end do
    if ( ma%has_bdybottom ) then
      do k = 1 , kz
        do j = jce1 , jce2
          vx(j,ice1,k) = v(j,ide1,k)
          vx(j,ici1,k) = 0.5_rkx * (v(j,ide1,k)+v(j,idi1,k))
        end do
      end do
    end if
    if ( ma%has_bdytop ) then
      do k = 1 , kz
        do j = jce1 , jce2
          vx(j,ice2,k) = 0.5_rkx*(v(j,ide2,k) + v(j,idi2,k))
        end do
      end do
    end if
  end subroutine uvstagtox

end module mod_moloch

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
