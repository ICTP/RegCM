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

module mod_advection
  !
  ! Horizontal and vertical advection.
  !
  use mod_regcm_types
  use mod_dynparam
  use mod_runparams
  use mod_memutil
  use mod_mpmessage
  use mod_constants
  use mod_service

  implicit none

  private

  public :: init_advection, hadv, vadv, start_advect

  real(rkx) :: ul

  interface hadv
    module procedure hadvuv
    module procedure hadvt
    module procedure hadvqv
    module procedure hadv3d
    module procedure hadvqx
    module procedure hadvtr
  end interface hadv

  interface vadv
    module procedure vadvuv
    module procedure vadv3d
    module procedure vadvqv
    module procedure vadv4d
  end interface vadv

  real(rkx), pointer, contiguous, dimension(:,:,:) :: ua => null( )   ! U/m * ps
  real(rkx), pointer, contiguous, dimension(:,:,:) :: va  => null( )  ! V/m * ps
  ! Surface pressure
  real(rkx), pointer, contiguous, dimension(:,:) :: ps  => null( )
  ! Surface pressure dot points
  real(rkx), pointer, contiguous, dimension(:,:) :: pd  => null( )
  ! 1/(mapfx**2 * 4 * dx)
  real(rkx), pointer, contiguous, dimension(:,:) :: xmapf  => null( )
  ! 1/(mapfd**2 * 16 * dx)
  real(rkx), pointer, contiguous, dimension(:,:) :: dmapf  => null( )
  ! Sigma Vertical Velocity
  real(rkx), pointer, contiguous, dimension(:,:,:) :: svv  => null( )
  ! Pressure full sigma levels
  real(rkx), pointer, contiguous, dimension(:,:,:) :: pfs  => null( )
  ! Pressure half sigma levels
  real(rkx), pointer, contiguous, dimension(:,:,:) :: phs  => null( )
  ! Mass divergence
  real(rkx), pointer, contiguous, dimension(:,:,:) :: divx  => null( )
  ! Top of PBL
  integer(ik4), pointer, contiguous, dimension(:,:) :: kpb  => null( )

  ! working space used to store the interlated values in vadv.

  real(rkx), pointer, contiguous, dimension(:) :: dds => null( )
  real(rkx), pointer, contiguous, dimension(:) :: xds => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: fg => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uavg1 => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uavg2 => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: vavg1 => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: vavg2 => null( )

  contains

    subroutine init_advection
      use mod_atm_interface, only : mddom, sfs, atms, atmx
      use mod_atm_interface, only : mdv, qdot, kpbl
      implicit none
      integer(ik4) :: k
      call assignpnt(atmx%umc,ua)
      call assignpnt(atmx%vmc,va)
      call assignpnt(sfs%psa,ps)
      call assignpnt(sfs%psdota,pd)
      call assignpnt(mddom%xmsf,xmapf)
      call assignpnt(mddom%dmsf,dmapf)
      call assignpnt(atms%pf3d,pfs)
      call assignpnt(atms%pb3d,phs)
      call assignpnt(mdv%cr,divx)
      call assignpnt(qdot,svv)
      call assignpnt(kpbl,kpb)
      call getmem1d(dds,1,kzp1,'mod_advection:dds')
      call getmem1d(xds,1,kz,'mod_advection:xds')
      call getmem3d(fg,jci1,jci2,ici1,ici2,1,kz,'advection:fg')
      call getmem3d(uavg1,jci1,jci2,ici1,ici2,1,kz,'advection:uavg1')
      call getmem3d(uavg2,jci1,jci2,ici1,ici2,1,kz,'advection:uavg2')
      call getmem3d(vavg1,jci1,jci2,ici1,ici2,1,kz,'advection:vavg1')
      call getmem3d(vavg2,jci1,jci2,ici1,ici2,1,kz,'advection:vavg2')
      xds(:) =  d_one / dsigma(:)
      dds(1) = d_zero
      dds(kzp1) = d_zero
      do k = 2, kz
        dds(k) = d_one / (dsigma(k) + dsigma(k-1))
      end do
      ul = uoffc * d_half * dt/dx
    end subroutine init_advection
    !
    ! Pre-compute
    !
    subroutine start_advect
      implicit none
      integer(ik4) :: i, j, k
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        uavg1(j,i,k) = ua(j,i+1,k)   + ua(j,i,k)
        uavg2(j,i,k) = ua(j+1,i+1,k) + ua(j+1,i,k)
        vavg1(j,i,k) = va(j+1,i,k)   + va(j,i,k)
        vavg2(j,i,k) = va(j+1,i+1,k) + va(j,i+1,k)
      end do
    end subroutine start_advect
    !
    ! UV advection
    !
    subroutine hadvuv(uten,vten,u,v)
      implicit none
      real(rkx), pointer, contiguous, intent (in), dimension(:,:,:) :: u, v
      real(rkx), pointer, contiguous, intent (inout), dimension(:,:,:) :: uten, vten

      real(rkx) :: ucmona, ucmonb, ucmonc, vcmona, vcmonb, vcmonc
      real(rkx) :: divd, diag, ff1, ff2, ff3, ff4
      integer(ik4) :: i, j, k
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'hadvuv'
      integer(ik4), save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      !
      ! ua, va : are p*u/m and p*v/m
      ! dmapf  : is 1/(mapfd**2 * 16dx)
      !
      if ( .not. upstream_mode ) then
        if ( idynamic == 1 ) then
          do k = 1, kz
            do i = idi1, idi2
              do j = jdi1, jdi2
                ucmona = ua(j,i+1,k)   + d_two*ua(j,i,k)   + ua(j,i-1,k)
                ucmonb = ua(j+1,i+1,k) + d_two*ua(j+1,i,k) + ua(j+1,i-1,k)
                ucmonc = ua(j-1,i+1,k) + d_two*ua(j-1,i,k) + ua(j-1,i-1,k)
                vcmona = va(j+1,i,k)   + d_two*va(j,i,k)   + va(j-1,i,k)
                vcmonb = va(j+1,i+1,k) + d_two*va(j,i+1,k) + va(j-1,i+1,k)
                vcmonc = va(j+1,i-1,k) + d_two*va(j,i-1,k) + va(j-1,i-1,k)
                ucmonb = ucmona + ucmonb
                ucmonc = ucmonc + ucmona
                vcmonb = vcmona + vcmonb
                vcmonc = vcmonc + vcmona
                uten(j,i,k) = uten(j,i,k) - dmapf(j,i) * &
                            ((u(j+1,i,k)+u(j,  i,k)) * ucmonb - &
                             (u(j,  i,k)+u(j-1,i,k)) * ucmonc + &
                             (u(j,i+1,k)+u(j,  i,k)) * vcmonb - &
                             (u(j,  i,k)+u(j,i-1,k)) * vcmonc)
                vten(j,i,k) = vten(j,i,k) - dmapf(j,i) * &
                            ((v(j+1,i,k)+v(j,  i,k)) * ucmonb - &
                             (v(j,  i,k)+v(j-1,i,k)) * ucmonc + &
                             (v(j,i+1,k)+v(j,  i,k)) * vcmonb - &
                             (v(j,  i,k)+v(j,i-1,k)) * vcmonc)
              end do
            end do
          end do
        else
          do k = 1, kz
            do i = idi1, idi2
              do j = jdi1, jdi2
                divd = d_rfour * ( divx(j,i,k)   + divx(j,i-1,k)   + &
                                   divx(j-1,i,k) + divx(j-1,i-1,k) )
                ucmona = ua(j,i+1,k)   + d_two*ua(j,i,k)   + ua(j,i-1,k)
                ucmonb = ua(j+1,i+1,k) + d_two*ua(j+1,i,k) + ua(j+1,i-1,k)
                ucmonc = ua(j-1,i+1,k) + d_two*ua(j-1,i,k) + ua(j-1,i-1,k)
                vcmona = va(j+1,i,k)   + d_two*va(j,i,k)   + va(j-1,i,k)
                vcmonb = va(j+1,i+1,k) + d_two*va(j,i+1,k) + va(j-1,i+1,k)
                vcmonc = va(j+1,i-1,k) + d_two*va(j,i-1,k) + va(j-1,i-1,k)
                diag = divd - dmapf(j,i) * &
                       ( ( ucmonb - ucmonc ) + ( vcmonb - vcmonc ) )
                ucmonb = ucmona + ucmonb
                ucmonc = ucmonc + ucmona
                vcmonb = vcmona + vcmonb
                vcmonc = vcmonc + vcmona
                uten(j,i,k) = uten(j,i,k) + u(j,i,k) * diag - dmapf(j,i) * &
                            (u(j+1,i,k) * ucmonb - u(j-1,i,k) * ucmonc + &
                             u(j,i+1,k) * vcmonb - u(j,i-1,k) * vcmonc)
                vten(j,i,k) = vten(j,i,k) + v(j,i,k) * diag - dmapf(j,i) * &
                            (v(j+1,i,k) * ucmonb - v(j-1,i,k) * ucmonc + &
                             v(j,i+1,k) * vcmonb - v(j,i-1,k) * vcmonc)
              end do
            end do
          end do
        end if
#ifdef DEBUG
        call time_end(subroutine_name,idindx)
#endif
        return
      end if

      if ( idynamic == 1 ) then
        do k = 1, kz
          do i = idi1, idi2
            do j = jdi1, jdi2
              ucmona = ua(j,i+1,k)   + d_two*ua(j,i,k)   + ua(j,i-1,k)
              ucmonb = ua(j+1,i+1,k) + d_two*ua(j+1,i,k) + ua(j+1,i-1,k)
              ucmonc = ua(j-1,i+1,k) + d_two*ua(j-1,i,k) + ua(j-1,i-1,k)
              vcmona = va(j+1,i,k)   + d_two*va(j,i,k)   + va(j-1,i,k)
              vcmonb = va(j+1,i+1,k) + d_two*va(j,i+1,k) + va(j-1,i+1,k)
              vcmonc = va(j+1,i-1,k) + d_two*va(j,i-1,k) + va(j-1,i-1,k)
              ff1 = ul*(u(j+1,i,k)+u(j,i,k))
              ff2 = ul*(u(j-1,i,k)+u(j,i,k))
              ff3 = ul*(v(j,i+1,k)+v(j,i,k))
              ff4 = ul*(v(j,i-1,k)+v(j,i,k))
              ucmonb = (d_one+ff1)*ucmona + (d_one-ff1)*ucmonb
              ucmonc = (d_one+ff2)*ucmonc + (d_one-ff2)*ucmona
              vcmonb = (d_one+ff3)*vcmona + (d_one-ff3)*vcmonb
              vcmonc = (d_one+ff4)*vcmonc + (d_one-ff4)*vcmona
              uten(j,i,k) = uten(j,i,k) - dmapf(j,i) * &
                          ((u(j+1,i,k)+u(j,  i,k)) * ucmonb - &
                           (u(j,  i,k)+u(j-1,i,k)) * ucmonc + &
                           (u(j,i+1,k)+u(j,  i,k)) * vcmonb - &
                           (u(j,  i,k)+u(j,i-1,k)) * vcmonc)
              vten(j,i,k) = vten(j,i,k) - dmapf(j,i) * &
                          ((v(j+1,i,k)+v(j,  i,k)) * ucmonb - &
                           (v(j,  i,k)+v(j-1,i,k)) * ucmonc + &
                           (v(j,i+1,k)+v(j,  i,k)) * vcmonb - &
                           (v(j,  i,k)+v(j,i-1,k)) * vcmonc)
            end do
          end do
        end do
      else
        do k = 1, kz
          do i = idi1, idi2
            do j = jdi1, jdi2
              divd = d_rfour * ( divx(j,i,k)   + divx(j,i-1,k)   + &
                                 divx(j-1,i,k) + divx(j-1,i-1,k) )
              ucmona = ua(j,i+1,k)   + d_two*ua(j,i,k)   + ua(j,i-1,k)
              ucmonb = ua(j+1,i+1,k) + d_two*ua(j+1,i,k) + ua(j+1,i-1,k)
              ucmonc = ua(j-1,i+1,k) + d_two*ua(j-1,i,k) + ua(j-1,i-1,k)
              vcmona = va(j+1,i,k)   + d_two*va(j,i,k)   + va(j-1,i,k)
              vcmonb = va(j+1,i+1,k) + d_two*va(j,i+1,k) + va(j-1,i+1,k)
              vcmonc = va(j+1,i-1,k) + d_two*va(j,i-1,k) + va(j-1,i-1,k)
              diag = divd - dmapf(j,i) * &
                      ( ( ucmonb - ucmonc ) + ( vcmonb - vcmonc ) )
              ff1 = ul*(u(j+1,i,k)+u(j,i,k))
              ff2 = ul*(u(j-1,i,k)+u(j,i,k))
              ff3 = ul*(v(j,i+1,k)+v(j,i,k))
              ff4 = ul*(v(j,i-1,k)+v(j,i,k))
              ucmonb = (d_one+ff1)*ucmona + (d_one-ff1)*ucmonb
              ucmonc = (d_one+ff2)*ucmonc + (d_one-ff2)*ucmona
              vcmonb = (d_one+ff3)*vcmona + (d_one-ff3)*vcmonb
              vcmonc = (d_one+ff4)*vcmonc + (d_one-ff4)*vcmona
              uten(j,i,k) = uten(j,i,k) + u(j,i,k) * diag - dmapf(j,i) * &
                          (u(j+1,i,k) * ucmonb - u(j-1,i,k) * ucmonc + &
                           u(j,i+1,k) * vcmonb - u(j,i-1,k) * vcmonc)
              vten(j,i,k) = vten(j,i,k) + v(j,i,k) * diag - dmapf(j,i) * &
                          (v(j+1,i,k) * ucmonb - v(j-1,i,k) * ucmonc + &
                           v(j,i+1,k) * vcmonb - v(j,i-1,k) * vcmonc)
            end do
          end do
        end do
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine hadvuv

    subroutine vadvuv(uten,vten,ud,vd)
      implicit none
      real(rkx), pointer, contiguous, intent (in), dimension(:,:,:) :: ud, vd
      real(rkx), pointer, contiguous, intent (inout), dimension(:,:,:) :: uten, vten

      real(rkx) :: qq, uu, vv
      integer(ik4) :: i, j, k
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'vadvuv'
      integer(ik4), save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      !
      ! vertical advection terms : interpolate ua or va to full sigma levels
      !
      do k = 2, kz
        do i = idi1, idi2
          do j = jdi1, jdi2
            qq = d_rfour * (svv(j,i,k)   + svv(j,i-1,k) + &
                            svv(j-1,i,k) + svv(j-1,i-1,k))
            uu = qq * (twt(k,1)*ud(j,i,k) + twt(k,2)*ud(j,i,k-1))
            vv = qq * (twt(k,1)*vd(j,i,k) + twt(k,2)*vd(j,i,k-1))
            uten(j,i,k-1) = uten(j,i,k-1) - uu*xds(k-1)
            uten(j,i,k)   = uten(j,i,k)   + uu*xds(k)
            vten(j,i,k-1) = vten(j,i,k-1) - vv*xds(k-1)
            vten(j,i,k)   = vten(j,i,k)   + vv*xds(k)
          end do
        end do
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine vadvuv
    !
    !  HADV
    !     This subroutines computes the horizontal flux-divergence terms.
    !     second-order difference is used.
    !     ften   : is the tendency for variable 'f'.
    !     f      : is p*f.
    !
    subroutine hadvt(ften,f)
      implicit none
      real(rkx), pointer, contiguous, intent (in), dimension(:,:,:) :: f
      real(rkx), pointer, contiguous, intent (inout), dimension(:,:,:) :: ften
      integer(ik4) :: i, j, k
      real(rkx) :: f1, f2, fx1, fx2, fy1, fy2
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'hadvt'
      integer(ik4), save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      if ( .not. upstream_mode ) then
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              fx1 = f(j-1,i,k) + f(j,i,k)
              fx2 = f(j,i,k)   + f(j+1,i,k)
              fy1 = f(j,i-1,k) + f(j,i,k)
              fy2 = f(j,i,k)   + f(j,i+1,k)
              fg(j,i,k) = - xmapf(j,i) * &
                     (uavg2(j,i,k)*fx2-uavg1(j,i,k)*fx1 + &
                      vavg2(j,i,k)*fy2-vavg1(j,i,k)*fy1)
            end do
          end do
        end do
      else
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              f1 = d_half*ul*(uavg2(j,i,k)+uavg1(j,i,k))/ps(j,i)
              f2 = d_half*ul*(vavg2(j,i,k)+vavg1(j,i,k))/ps(j,i)
              fx1 = (d_one+f1)*f(j-1,i,k) + (d_one-f1)*f(j,i,k)
              fx2 = (d_one+f1)*f(j,i,k)   + (d_one-f1)*f(j+1,i,k)
              fy1 = (d_one+f2)*f(j,i-1,k) + (d_one-f2)*f(j,i,k)
              fy2 = (d_one+f2)*f(j,i,k)   + (d_one-f2)*f(j,i+1,k)
              fg(j,i,k) = - xmapf(j,i) * &
                     (uavg2(j,i,k)*fx2 - uavg1(j,i,k)*fx1 + &
                      vavg2(j,i,k)*fy2 - vavg1(j,i,k)*fy1)
            end do
          end do
        end do
      end if
      !
      ! Instability correction
      !
      ! Local extrema exceeding a certain threshold
      ! must not grow further due to advection
      !
      if ( stability_enhance ) then
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( abs(f(j,i+1,k) + f(j,i-1,k) - &
                       d_two*f(j,i,k))/ps(j,i) > t_extrema ) then
                if ( (f(j,i,k) > f(j,i+1,k)) .and. &
                     (f(j,i,k) > f(j,i-1,k)) ) then
                  fg(j,i,k) = min(fg(j,i,k),d_zero)
                else if ( (f(j,i,k) < f(j,i+1,k)) .and. &
                          (f(j,i,k) < f(j,i-1,k)) ) then
                  fg(j,i,k) = max(fg(j,i,k),d_zero)
                end if
              end if
              if ( abs(f(j+1,i,k) + f(j-1,i,k) - &
                       d_two*f(j,i,k))/ps(j,i) > t_extrema ) then
                if ( (f(j,i,k) > f(j+1,i,k)) .and. &
                     (f(j,i,k) > f(j-1,i,k)) ) then
                  fg(j,i,k) = min(fg(j,i,k),d_zero)
                else if ( (f(j,i,k) < f(j+1,i,k)) .and. &
                          (f(j,i,k) < f(j-1,i,k)) ) then
                  fg(j,i,k) = max(fg(j,i,k),d_zero)
                end if
              end if
            end do
          end do
        end do
      end if
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        ften(j,i,k) = ften(j,i,k) + fg(j,i,k)
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine hadvt

    subroutine hadv3d(ften,f,ind)
      implicit none
      integer(ik4), intent (in) :: ind
      real(rkx), pointer, contiguous, intent (in), dimension(:,:,:) :: f
      real(rkx), pointer, contiguous, intent (inout), dimension(:,:,:) :: ften
      integer(ik4) :: i, j, k
      real(rkx) :: uaz1, uaz2, vaz1, vaz2
      real(rkx) :: f1, f2, fx1, fx2, fy1, fy2
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'hadv3d'
      integer(ik4), save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif

      if ( .not. upstream_mode ) then
        if ( ind == 0 ) then
          !
          ! for cross point variables on half sigma levels
          !
          do k = 1, kz
            do i = ici1, ici2
              do j = jci1, jci2
                fx1 = f(j-1,i,k) + f(j,i,k)
                fx2 = f(j,i,k)   + f(j+1,i,k)
                fy1 = f(j,i-1,k) + f(j,i,k)
                fy2 = f(j,i,k)   + f(j,i+1,k)
                ften(j,i,k) = ften(j,i,k) - xmapf(j,i) * &
                              (uavg2(j,i,k)*fx2 - uavg1(j,i,k)*fx1 + &
                               vavg2(j,i,k)*fy2 - vavg1(j,i,k)*fy1)
              end do
            end do
          end do
        else if ( ind == 1 ) then
          !
          ! Interpolate the winds to the full sigma levels
          ! while the advection term is calculated
          !
          do k = 2, kz
            do i = ici1, ici2
              do j = jci1, jci2
                uaz1 = ( twt(k,1) * uavg1(j,i,k) + &
                         twt(k,2) * uavg1(j,i,k-1) )
                uaz2 = ( twt(k,1) * uavg2(j,i,k) + &
                         twt(k,2) * uavg2(j,i,k-1) )
                vaz1 = ( twt(k,1) * vavg1(j,i,k) + &
                         twt(k,2) * vavg1(j,i,k-1) )
                vaz2 = ( twt(k,1) * vavg2(j,i,k) + &
                         twt(k,2) * vavg2(j,i,k-1) )
                fx1 = f(j-1,i,k) + f(j,i,k)
                fx2 = f(j,i,k)   + f(j+1,i,k)
                fy1 = f(j,i-1,k) + f(j,i,k)
                fy2 = f(j,i,k)   + f(j,i+1,k)
                ften(j,i,k) = ften(j,i,k) - xmapf(j,i) * &
                              (uaz2*fx2-uaz1*fx1 + vaz2*fy2-vaz1*fy1)
              end do
            end do
          end do
        else
          call fatal(__FILE__,__LINE__, &
                     'Unsupported advection scheme')
        end if
#ifdef DEBUG
        call time_end(subroutine_name,idindx)
#endif
        return
      end if

      if ( ind == 0 ) then
        !
        ! for cross point variables on half sigma levels
        !
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              f1 = d_half*ul*(uavg2(j,i,k)+uavg1(j,i,k))/ps(j,i)
              f2 = d_half*ul*(vavg2(j,i,k)+vavg1(j,i,k))/ps(j,i)
              fx1 = (d_one+f1)*f(j-1,i,k) + (d_one-f1)*f(j,i,k)
              fx2 = (d_one+f1)*f(j,i,k)   + (d_one-f1)*f(j+1,i,k)
              fy1 = (d_one+f2)*f(j,i-1,k) + (d_one-f2)*f(j,i,k)
              fy2 = (d_one+f2)*f(j,i,k)   + (d_one-f2)*f(j,i+1,k)
              ften(j,i,k) = ften(j,i,k) - xmapf(j,i) * &
                        (uavg2(j,i,k)*fx2 - uavg1(j,i,k)*fx1 + &
                         vavg2(j,i,k)*fy2 - vavg1(j,i,k)*fy1)
            end do
          end do
        end do
      else if ( ind == 1 ) then
        !
        ! Interpolate the winds to the full sigma levels
        ! while the advection term is calculated
        !
        do k = 2, kz
          do i = ici1, ici2
            do j = jci1, jci2
              uaz1 = ( twt(k,1) * uavg1(j,i,k) + &
                       twt(k,2) * uavg1(j,i,k-1) )
              uaz2 = ( twt(k,1) * uavg2(j,i,k) + &
                       twt(k,2) * uavg2(j,i,k-1) )
              vaz1 = ( twt(k,1) * vavg1(j,i,k) + &
                       twt(k,2) * vavg1(j,i,k-1) )
              vaz2 = ( twt(k,1) * vavg2(j,i,k) + &
                       twt(k,2) * vavg2(j,i,k-1) )
              f1 = d_half*ul*(uavg2(j,i,k)+uavg1(j,i,k))/ps(j,i)
              f2 = d_half*ul*(vavg2(j,i,k)+vavg1(j,i,k))/ps(j,i)
              fx1 = (d_one+f1)*f(j-1,i,k) + (d_one-f1)*f(j,i,k)
              fx2 = (d_one+f1)*f(j,i,k)   + (d_one-f1)*f(j+1,i,k)
              fy1 = (d_one+f2)*f(j,i-1,k) + (d_one-f2)*f(j,i,k)
              fy2 = (d_one+f2)*f(j,i,k)   + (d_one-f2)*f(j,i+1,k)
              ften(j,i,k) = ften(j,i,k) - xmapf(j,i) * &
                        (uaz2*fx2-uaz1*fx1 + vaz2*fy2-vaz1*fy1)
            end do
          end do
        end do
      else
        call fatal(__FILE__,__LINE__, &
                   'Unsupported advection scheme')
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine hadv3d

    subroutine hadvqv(ften,f,iv)
      implicit none
      integer(ik4), intent(in) :: iv
      real(rkx), pointer, contiguous, intent (in), dimension(:,:,:,:) :: f
      real(rkx), pointer, contiguous, intent (inout), dimension(:,:,:,:) :: ften
      real(rkx) :: f1, f2, fx1, fx2, fy1, fy2
      integer(ik4) :: i, j, k
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'hadvqv'
      integer(ik4), save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      !
      ! for qv:
      !
      if ( .not. upstream_mode ) then
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              fx1 = f(j-1,i,k,iv) + f(j,i,k,iv)
              fx2 = f(j,i,k,iv)   + f(j+1,i,k,iv)
              fy1 = f(j,i-1,k,iv) + f(j,i,k,iv)
              fy2 = f(j,i,k,iv)   + f(j,i+1,k,iv)
              fg(j,i,k) = -xmapf(j,i) * &
                     (uavg2(j,i,k)*fx2 - uavg1(j,i,k)*fx1 + &
                      vavg2(j,i,k)*fy2 - vavg1(j,i,k)*fy1)
            end do
          end do
        end do
      else
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              f1 = d_half*ul*(uavg2(j,i,k)+uavg1(j,i,k))/ps(j,i)
              f2 = d_half*ul*(vavg2(j,i,k)+vavg1(j,i,k))/ps(j,i)
              fx1 = (d_one+f1)*f(j-1,i,k,iv) + (d_one-f1)*f(j,i,k,iv)
              fx2 = (d_one+f1)*f(j,i,k,iv)   + (d_one-f1)*f(j+1,i,k,iv)
              fy1 = (d_one+f2)*f(j,i-1,k,iv) + (d_one-f2)*f(j,i,k,iv)
              fy2 = (d_one+f2)*f(j,i,k,iv)   + (d_one-f2)*f(j,i+1,k,iv)
              fg(j,i,k) = -xmapf(j,i) * &
                        (uavg2(j,i,k)*fx2 - uavg1(j,i,k)*fx1 + &
                         vavg2(j,i,k)*fy2 - vavg1(j,i,k)*fy1)
            end do
          end do
        end do
      end if
      !
      ! Instability correction
      !
      ! Local extrema exceeding a certain realtive threshold
      ! must not grow further due to advection
      !
      if ( stability_enhance ) then
        do k = 1, kz
          do i = ici1, ici2
            do j = jci1, jci2
              if ( abs(f(j,i+1,k,iv) + f(j,i-1,k,iv) - d_two*f(j,i,k,iv)) / &
                      max(f(j,i,k,iv),dlowval) > q_rel_extrema ) then
                if ( (f(j,i,k,iv) > f(j,i+1,k,iv)) .and. &
                     (f(j,i,k,iv) > f(j,i-1,k,iv)) ) then
                  fg(j,i,k) = min(fg(j,i,k),d_zero)
                else if ( (f(j,i,k,iv) < f(j,i+1,k,iv)) .and. &
                          (f(j,i,k,iv) < f(j,i-1,k,iv)) ) then
                  fg(j,i,k) = max(fg(j,i,k),d_zero)
                end if
              end if
              if ( abs(f(j+1,i,k,iv) + f(j-1,i,k,iv) - d_two*f(j,i,k,iv)) / &
                      max(f(j,i,k,iv),dlowval) > q_rel_extrema ) then
                if ( (f(j,i,k,iv) > f(j+1,i,k,iv)) .and. &
                     (f(j,i,k,iv) > f(j-1,i,k,iv)) ) then
                  fg(j,i,k) = min(fg(j,i,k),d_zero)
                else if ( (f(j,i,k,iv) < f(j+1,i,k,iv)) .and. &
                          (f(j,i,k,iv) < f(j-1,i,k,iv)) ) then
                  fg(j,i,k) = max(fg(j,i,k),d_zero)
                end if
              end if
            end do
          end do
        end do
      end if
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        ften(j,i,k,iv) = ften(j,i,k,iv) + fg(j,i,k)
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine hadvqv
    !
    ! Up-wind values are used for non-hydrostatic
    !
    subroutine hadvqx(ften,f,n1,n2)
      implicit none
      integer(ik4), intent(in) :: n1, n2
      real(rkx), pointer, contiguous, intent (in), dimension(:,:,:,:) :: f
      real(rkx), pointer, contiguous, intent (inout), dimension(:,:,:,:) :: ften

      integer(ik4) :: i, j, k, n
      real(rkx) :: f1, f2, fx1, fx2, fy1, fy2
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'hadvqx'
      integer(ik4), save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      !
      ! for qx different from qv
      !
      do n = n1, n2
        if ( .not. upstream_mode ) then
          do k = 1, kz
            do i = ici1, ici2
              do j = jci1, jci2
                fx1 = f(j-1,i,k,n) + f(j,i,k,n)
                fx2 = f(j,i,k,n)   + f(j+1,i,k,n)
                fy1 = f(j,i-1,k,n) + f(j,i,k,n)
                fy2 = f(j,i,k,n)   + f(j,i+1,k,n)
                fg(j,i,k) = - xmapf(j,i) * &
                     (uavg2(j,i,k)*fx2 - uavg1(j,i,k)*fx1 + &
                      vavg2(j,i,k)*fy2 - vavg1(j,i,k)*fy1)
              end do
            end do
          end do
        else
          do k = 1, kz
            do i = ici1, ici2
              do j = jci1, jci2
                f1 = d_half*ul*(uavg2(j,i,k)+uavg1(j,i,k))/ps(j,i)
                f2 = d_half*ul*(vavg2(j,i,k)+vavg1(j,i,k))/ps(j,i)
                fx1 = (d_one+f1)*f(j-1,i,k,n) + (d_one-f1)*f(j,i,k,n)
                fx2 = (d_one+f1)*f(j,i,k,n)   + (d_one-f1)*f(j+1,i,k,n)
                fy1 = (d_one+f2)*f(j,i-1,k,n) + (d_one-f2)*f(j,i,k,n)
                fy2 = (d_one+f2)*f(j,i,k,n)   + (d_one-f2)*f(j,i+1,k,n)
                fg(j,i,k) = - xmapf(j,i) * &
                           (uavg2(j,i,k)*fx2 - uavg1(j,i,k)*fx1 + &
                            vavg2(j,i,k)*fy2 - vavg1(j,i,k)*fy1)
              end do
            end do
          end do
        end if
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          ften(j,i,k,n) = ften(j,i,k,n) + fg(j,i,k)
        end do
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine hadvqx
    !
    ! Up-wind values are used
    !
    subroutine hadvtr(ften,f)
      implicit none
      real(rkx), pointer, contiguous, intent (in), dimension(:,:,:,:) :: f
      real(rkx), pointer, contiguous, intent (inout), dimension(:,:,:,:) :: ften

      integer(ik4) :: i, j, k, n
      real(rkx) :: f1, f2, fx1, fx2, fy1, fy2
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'hadvtr'
      integer(ik4), save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      !
      ! for tracers
      !
      do n = 1, ntr
        if ( .not. upstream_mode ) then
          do k = 1, kz
            do i = ici1, ici2
              do j = jci1, jci2
                fx1 = f(j-1,i,k,n)+f(j,i,k,n)
                fx2 = f(j,i,k,n)  +f(j+1,i,k,n)
                fy1 = f(j,i-1,k,n)+f(j,i,k,n)
                fy2 = f(j,i,k,n)  +f(j,i+1,k,n)
                fg(j,i,k) = - xmapf(j,i) * &
                        (uavg2(j,i,k)*fx2 - uavg1(j,i,k)*fx1 + &
                         vavg2(j,i,k)*fy2 - vavg1(j,i,k)*fy1)
              end do
            end do
          end do
        else
          do k = 1, kz
            do i = ici1, ici2
              do j = jci1, jci2
                f1 = d_half*ul*(uavg2(j,i,k)+uavg1(j,i,k))/ps(j,i)
                f2 = d_half*ul*(vavg2(j,i,k)+vavg1(j,i,k))/ps(j,i)
                fx1 = (d_one+f1)*f(j-1,i,k,n) + (d_one-f1)*f(j,i,k,n)
                fx2 = (d_one+f1)*f(j,i,k,n)   + (d_one-f1)*f(j+1,i,k,n)
                fy1 = (d_one+f2)*f(j,i-1,k,n) + (d_one-f2)*f(j,i,k,n)
                fy2 = (d_one+f2)*f(j,i,k,n)   + (d_one-f2)*f(j,i+1,k,n)
                fg(j,i,k) = - xmapf(j,i) * &
                    (uavg2(j,i,k)*fx2 - uavg1(j,i,k)*fx1 + &
                     vavg2(j,i,k)*fy2 - vavg1(j,i,k)*fy1)
              end do
            end do
          end do
        end if
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          ften(j,i,k,n) = ften(j,i,k,n) + fg(j,i,k)
        end do
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine hadvtr
    !
    ! VADV
    !     This subroutine computes the vertical flux-divergence terms.
    !     ften   : is the tendency of variable 'f'.
    !     f      : is p*f.
    !     ind = 0 : for pp, w
    !           1 : for t.
    !           2 : Use pbl information
    !
    subroutine vadv3d(ften,f,nk,ind)
      implicit none
      integer(ik4), intent(in) :: ind, nk
      real(rkx), pointer, contiguous, intent (in), dimension(:,:,:) :: f
      real(rkx), pointer, contiguous, intent (inout), dimension(:,:,:) :: ften

      real(rkx) :: rdphf, rdplf, fx, qq
      real(rkx), dimension(jci1:jci2,ici1:ici2,kz) :: dotqdot
      integer(ik4) :: i, j, k
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'vadv3d'
      integer(ik4), save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      if ( ind == 0 ) then
        if ( nk == kz ) then
          do k = 2, kz
            do i = ici1, ici2
              do j = jci1, jci2
                fx = svv(j,i,k) * (twt(k,1)*f(j,i,k)+twt(k,2)*f(j,i,k-1))
                ften(j,i,k-1) = ften(j,i,k-1) - fx*xds(k-1)
                ften(j,i,k)   = ften(j,i,k)   + fx*xds(k)
              end do
            end do
          end do
        else
          do k = 1, kz
            do i = ici1, ici2
              do j = jci1, jci2
                qq = d_half * (svv(j,i,k) + svv(j,i,k+1))
                fx = qq * ((f(j,i,k) + f(j,i,k+1)))
                ften(j,i,k+1) = ften(j,i,k+1) + fx*dds(k+1)
                ften(j,i,k)   = ften(j,i,k)   - fx*dds(k)
              end do
            end do
          end do
        end if
      else if ( ind == 1 ) then
        !
        ! vertical advection terms : interpolate to full sigma levels
        !
        if ( idynamic == 1 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2 )
            dotqdot(j,i,1) = d_zero
          end do
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:nk )
            dotqdot(j,i,k) = svv(j,i,k) * &
                    (twt(k,1)*f(j,i,k)  *(pfs(j,i,k)/phs(j,i,k))  **c287 + &
                     twt(k,2)*f(j,i,k-1)*(pfs(j,i,k)/phs(j,i,k-1))**c287)
          end do
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:nk )
            ften(j,i,k-1) = ften(j,i,k-1) - dotqdot(j,i,k)*xds(k-1)
            ften(j,i,k)   = ften(j,i,k)   + dotqdot(j,i,k)*xds(k)
          end do
        else
          do i = ici1, ici2
            do j = jci1, jci2
              rdphf = exp(-c287*log(phs(j,i,1)))
              dotqdot(j,i,1) = f(j,i,1)*rdphf
            end do
          end do
          do k = 2, kz
            do i = ici1, ici2
              do j = jci1, jci2
                rdphf = exp(-c287*log(phs(j,i,k)))
                rdplf = exp( c287*log(pfs(j,i,k)))
                dotqdot(j,i,k) = f(j,i,k)*rdphf
                fx = rdplf * svv(j,i,k) * ( twt(k,1) * dotqdot(j,i,k) + &
                                            twt(k,2) * dotqdot(j,i,k-1) )
                ften(j,i,k-1) = ften(j,i,k-1) - fx*xds(k-1)
                ften(j,i,k)   = ften(j,i,k)   + fx*xds(k)
              end do
            end do
          end do
        end if
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine vadv3d

    subroutine vadvqv(ften,f,n)
      implicit none
      integer(ik4), intent(in) :: n
      real(rkx), pointer, contiguous, intent (in), dimension(:,:,:,:) :: f
      real(rkx), pointer, contiguous, intent (inout), dimension(:,:,:,:) :: ften
      integer(ik4) :: i, j, k
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'vadvqv'
      integer(ik4), save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        fg(j,i,k) = d_zero
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
        if ( f(j,i,k,n)   > minqq * ps(j,i) .and. &
             f(j,i,k-1,n) > minqq * ps(j,i) ) then
          fg(j,i,k) = f(j,i,k,n)*(f(j,i,k-1,n)/f(j,i,k,n))**qcon(k)
        end if
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
        ften(j,i,k-1,n) = ften(j,i,k-1,n)-svv(j,i,k)*fg(j,i,k)*xds(k-1)
        ften(j,i,k,n)   = ften(j,i,k,n)  +svv(j,i,k)*fg(j,i,k)*xds(k)
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine vadvqv
    !
    ! VADV
    !     This subroutine computes the vertical flux-divergence terms.
    !     ften   : is the tendency of variable 'f'.
    !     f      : is p*f.
    !     ind = 1 : for qv
    !           2 : for chemical tracers
    !           3 : for hydometeors
    !           4 : use pbl information
    !
    subroutine vadv4d(ften,f,n1,n2,ind)
      implicit none
      integer(ik4), intent(in) :: ind, n1, n2
      real(rkx), pointer, contiguous, intent (in), dimension(:,:,:,:) :: f
      real(rkx), pointer, contiguous, intent (inout), dimension(:,:,:,:) :: ften
      real(rkx) :: slope
      integer(ik4) :: i, j, k, n
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'vadv4d'
      integer(ik4), save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      do n = n1, n2
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          fg(j,i,k) = d_zero
        end do
        if ( ind == 0 ) then
          do k = 2, kz
            do i = ici1, ici2
              do j = jci1, jci2
                if ( svv(j,i,k) > d_zero ) then
                  fg(j,i,k) = svv(j,i,k)*f(j,i,k-1,n)
                else
                  fg(j,i,k) = svv(j,i,k)*f(j,i,k,n)
                end if
              end do
            end do
          end do
        else if ( ind == 1 ) then
          do k = 2, kz
            do i = ici1, ici2
              do j = jci1, jci2
                if ( svv(j,i,k) > d_zero ) then
                  if ( f(j,i,k-1,n) > minqq * minqq * ps(j,i)  ) then
                    fg(j,i,k) = svv(j,i,k) * &
                        (twt(k,1)*f(j,i,k,n) + twt(k,2)*f(j,i,k-1,n))
                  else
                    fg(j,i,k) = d_zero
                  end if
                else
                  if ( f(j,i,k,n) > minqq * minqq  * ps(j,i) ) then
                    fg(j,i,k) = svv(j,i,k) * &
                        (twt(k,1)*f(j,i,k,n) + twt(k,2)*f(j,i,k-1,n))
                  else
                    fg(j,i,k) = d_zero
                  end if
                end if
              end do
            end do
          end do
        else if ( ind == 2 ) then
          do k = 2, kz
            do i = ici1, ici2
              do j = jci1, jci2
                if ( svv(j,i,k) > d_zero ) then
                  if ( f(j,i,k-1,n) > mintr * ps(j,i)  ) then
                    fg(j,i,k) = svv(j,i,k) * &
                        (twt(k,1)*f(j,i,k,n) + twt(k,2)*f(j,i,k-1,n))
                  else
                    fg(j,i,k) = d_zero
                  end if
                else
                  if ( f(j,i,k,n) > mintr * ps(j,i) ) then
                    fg(j,i,k) = svv(j,i,k) * &
                        (twt(k,1)*f(j,i,k,n) + twt(k,2)*f(j,i,k-1,n))
                  else
                    fg(j,i,k) = d_zero
                  end if
                end if
              end do
            end do
          end do
        else if ( ind == 3 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
            fg(j,i,k) = twt(k,1)*f(j,i,k,n) + twt(k,2)*f(j,i,k-1,n)
          end do
          do i = ici1, ici2
            do j = jci1, jci2
              if ( kpb(j,i) > kz ) then
                call fatal(__FILE__,__LINE__,'kpbl is greater than kz')
              end if
              if ( kpb(j,i) >= 4 ) then
                ! Calculate slope of scalar in layer above ambiguous layer
                k = kpb(j,i)-2
                if ( (f(j,i,k+1,n)-f(j,i,k,n)) > d_zero .and. &
                     (f(j,i,k,n)-f(j,i,k-1,n)) > d_zero ) then
                  slope = min((f(j,i,k+1,n)-f(j,i,k,n)) / &
                          (hsigma(k+1)-hsigma(k)), &
                          (f(j,i,k,n)-f(j,i,k-1,n))/(hsigma(k)-hsigma(k-1)))
                else if ( (f(j,i,k+1,n)-f(j,i,k,n)) < d_zero .and. &
                          (f(j,i,k,n)-f(j,i,k-1,n)) < d_zero ) then
                  slope = max((f(j,i,k+1,n)-f(j,i,k,n)) / &
                          (hsigma(k+1)-hsigma(k)), &
                          (f(j,i,k,n)-f(j,i,k-1,n))/(hsigma(k)-hsigma(k-1)))
                else
                  slope = d_zero
                end if
                ! Replace the values of scalar at top and bottom of ambiguous
                ! layer as long as inversion is actually in the ambiguous layer
                k = kpb(j,i)
                fg(j,i,k-1) = f(j,i,k-2,n) + slope*(sigma(k-1)-hsigma(k-2))
                if (abs(f(j,i,k-2,n) + &
                        slope*(hsigma(k-1)-hsigma(k-2))-f(j,i,k,n)) > &
                    abs(f(j,i,k-1,n)-f(j,i,k,n)) ) then
                  fg(j,i,k) = f(j,i,k,n)
                else
                  fg(j,i,k) = f(j,i,k-2,n) + slope*(sigma(k)-hsigma(k-2))
                end if
              end if
            end do
          end do
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
            fg(j,i,k) = fg(j,i,k) * svv(j,i,k)
          end do
        end if
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 2:kz )
          ften(j,i,k-1,n) = ften(j,i,k-1,n)-fg(j,i,k)*xds(k-1)
          ften(j,i,k,n)   = ften(j,i,k,n)  +fg(j,i,k)*xds(k)
        end do
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine vadv4d

end module mod_advection

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
