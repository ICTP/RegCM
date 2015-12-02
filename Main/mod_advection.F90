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

  public :: init_advection, hadv , vadv

  interface hadv
    module procedure hadv3d
    module procedure hadvuv
    module procedure hadv4d
  end interface hadv

  interface vadv
    module procedure vadv3d
    module procedure vadv4d
    module procedure vadvuv
  end interface vadv

  real(rk8) , pointer , dimension(:,:,:) :: ua   ! U wind * ps
  real(rk8) , pointer , dimension(:,:,:) :: va   ! V wind * ps
  real(rk8) , pointer , dimension(:,:) :: ps     ! Surface pressure
  real(rk8) , pointer , dimension(:,:) :: mapfx  ! Map factor Cross
  real(rk8) , pointer , dimension(:,:) :: mapfd  ! Map factor Dot
  real(rk8) , pointer , dimension(:,:) :: xmapf  ! 1/(mapfx**2*dx4)
  real(rk8) , pointer , dimension(:,:) :: dmapf  ! 1/(mapfd**2*dx16)
  real(rk8) , pointer , dimension(:,:,:) :: svv  ! Sigma Vertical Velocity
  real(rk8) , pointer , dimension(:,:,:) :: pfs  ! Pressure full sigma levels
  real(rk8) , pointer , dimension(:,:,:) :: phs  ! Pressure half sigma levels
  real(rk8) , pointer , dimension(:,:,:) :: divx ! Mass divergence
  integer(ik4) , pointer , dimension(:,:) :: kpb ! Top of PBL

  ! working space used to store the interlated values in vadv.

  real(rk8) , pointer , dimension(:) :: dds , xds
  real(rk8) , pointer , dimension(:,:,:) :: fg

  contains

    subroutine init_advection
      use mod_atm_interface , only : mddom , sfs , atms , atmc
      use mod_atm_interface , only : mdv , qdot , kpbl
      implicit none
      integer(ik4) :: k
      call assignpnt(atmc%u,ua)
      call assignpnt(atmc%v,va)
      call assignpnt(sfs%psa,ps)
      call assignpnt(mddom%msfx,mapfx)
      call assignpnt(mddom%msfd,mapfd)
      call assignpnt(mddom%xmsf,xmapf)
      call assignpnt(mddom%dmsf,dmapf)
      call assignpnt(atms%pf3d,pfs)
      call assignpnt(atms%pb3d,phs)
      call assignpnt(mdv%cr,divx)
      call assignpnt(qdot,svv)
      call assignpnt(kpbl,kpb)
      call getmem1d(dds,1,kzp1,'mod_advection:dds')
      call getmem1d(xds,1,kz,'mod_advection:xds')
      if ( ibltyp == 2 .and. iuwvadv == 1 ) then
        call getmem3d(fg,jci1,jci2,ici1,ici2,1,kz,'advection:fg')
      end if
      xds(:) =  d_one / dsigma(:)
      dds(1) = d_zero
      dds(kzp1) = d_zero
      do k = 2 , kz
        dds(k) = d_one / (dsigma(k) + dsigma(k-1))
      end do
    end subroutine init_advection
    !
    ! UV advection
    !
    subroutine hadvuv(uten,vten,u,v)
      implicit none
      real(rk8) , pointer , intent (in) , dimension(:,:,:) :: u , v
      real(rk8) , pointer , intent (inout), dimension(:,:,:) :: uten , vten

      real(rk8) :: ucmona , ucmonb , ucmonc , vcmona , vcmonb , vcmonc
      real(rk8) :: divd , diag
      integer(ik4) :: i , j , k
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'hadvuv'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      !
      ! ua, va : are p*u and p*v.
      ! msfd   : is the map scale factor at dot points.
      !
      if ( idynamic == 2 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
              divd = d_rfour * ( divx(j,i,k)   + divx(j,i-1,k)   + &
                                 divx(j-1,i,k) + divx(j-1,i-1,k) )
              ucmona = ua(j,i+1,k)   + d_two*ua(j,i,k)   + ua(j,i-1,k)
              vcmona = va(j+1,i,k)   + d_two*va(j,i,k)   + va(j-1,i,k)
              ucmonb = ua(j+1,i+1,k) + d_two*ua(j+1,i,k) + ua(j+1,i-1,k)
              vcmonb = va(j+1,i+1,k) + d_two*va(j,i+1,k) + va(j-1,i+1,k)
              ucmonc = ua(j-1,i+1,k) + d_two*ua(j-1,i,k) + ua(j-1,i-1,k)
              vcmonc = va(j+1,i-1,k) + d_two*va(j,i-1,k) + va(j-1,i-1,k)
              diag = divd - dmapf(j,i)*( ( ucmonb - ucmonc ) + &
                                         ( vcmonb - vcmonc ) )
              ucmonb = ucmonb + ucmona
              vcmonb = vcmonb + vcmona
              ucmonc = ucmonc + ucmona
              vcmonc = vcmonc + vcmona
              uten(j,i,k) = uten(j,i,k) + u(j,i,k) * diag -    &
                dmapf(j,i) * (u(j+1,i,k) * ucmonb - &
                              u(j-1,i,k) * ucmonc + &
                              u(j,i+1,k) * vcmonb - &
                              u(j,i-1,k) * vcmonc)
              vten(j,i,k) = vten(j,i,k) + v(j,i,k) * diag -    &
                dmapf(j,i) * (v(j+1,i,k) * ucmonb - &
                              v(j-1,i,k) * ucmonc + &
                              v(j,i+1,k) * vcmonb - &
                              v(j,i-1,k) * vcmonc)
            end do
          end do
        end do
      else
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
              ucmona = ua(j,i+1,k)+d_two*ua(j,i,k)+ua(j,i-1,k)
              vcmona = va(j+1,i,k)+d_two*va(j,i,k)+va(j-1,i,k)
              ucmonb = ua(j+1,i+1,k) + d_two*ua(j+1,i,k) + &
                       ua(j+1,i-1,k) + ucmona
              vcmonb = va(j+1,i+1,k) + d_two*va(j,i+1,k) + &
                       va(j-1,i+1,k) + vcmona
              ucmonc = ua(j-1,i+1,k) + d_two*ua(j-1,i,k) + &
                       ua(j-1,i-1,k) + ucmona
              vcmonc = va(j+1,i-1,k) + d_two*va(j,i-1,k) + &
                       va(j-1,i-1,k) + vcmona
              uten(j,i,k) = uten(j,i,k) - dmapf(j,i) *      &
                          ((u(j+1,i,k)+u(j,i,k))  *ucmonb - &
                           (u(j,i,k)  +u(j-1,i,k))*ucmonc + &
                           (u(j,i+1,k)+u(j,i,k))  *vcmonb - &
                           (u(j,i,k)  +u(j,i-1,k))*vcmonc)
              vten(j,i,k) = vten(j,i,k) - dmapf(j,i) *      &
                          ((v(j+1,i,k)+v(j,i,k))  *ucmonb - &
                           (v(j,i,k)  +v(j-1,i,k))*ucmonc + &
                           (v(j,i+1,k)+v(j,i,k))  *vcmonb - &
                           (v(j,i,k)  +v(j,i-1,k))*vcmonc)
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
      real(rk8) , pointer , intent (in) , dimension(:,:,:) :: ud , vd
      real(rk8) , pointer , intent (inout), dimension(:,:,:) :: uten , vten

      real(rk8) :: qq , uu , vv
      integer(ik4) :: i , j , k
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'vadvuv'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      !
      ! vertical advection terms : interpolate ua or va to full sigma levels
      !
      do k = 2 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            qq = d_rfour * (svv(j,i,k) + svv(j,i-1,k) + &
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
    !     ldot   : cross/dot variable flag
    !     ften   : is the tendency for variable 'f'.
    !     f      : is p*f.
    !     nk     : is the number of vertical levels to work (kz/kzp1)
    !
    subroutine hadv3d(ften,f,nk)
      implicit none
      integer(ik4) , intent (in) :: nk
      real(rk8) , pointer , intent (in) , dimension(:,:,:) :: f
      real(rk8) , pointer , intent (inout), dimension(:,:,:) :: ften

      integer(ik4) :: i , j , k
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'hadv3d'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      !
      ! for cross point variables on half sigma levels
      !
      if ( nk == kz ) then
        do k = 1 , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,k) = ften(j,i,k) - xmapf(j,i) *                  &
                  ((ua(j+1,i+1,k)+ua(j+1,i,k))*(f(j+1,i,k)+f(j,i,k)) -  &
                   (ua(j,i+1,k)+ua(j,i,k)) *   (f(j,i,k)+f(j-1,i,k)) +  &
                   (va(j+1,i+1,k)+va(j,i+1,k))*(f(j,i+1,k)+f(j,i,k)) -  &
                   (va(j+1,i,k)+va(j,i,k)) *   (f(j,i-1,k)+f(j,i,k)))
            end do
          end do
        end do
      else
        !
        ! Interpolate the winds to the full sigma levels
        ! while the advection term is calculated
        !
        do k = 2 , nk - 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,k) = ften(j,i,k) - xmapf(j,i)       *  &
                (((ua(j+1,i+1,k-1)+ua(j+1,i,k-1))*twt(k,2) +  &
                  (ua(j+1,i+1,k)  +ua(j+1,i,k))*twt(k,1)) *   &
                  (f(j+1,i,k)+f(j,i,k)) -                     &
                 ((ua(j,i+1,k-1)+ua(j,i,k-1))*twt(k,2) +      &
                  (ua(j,i+1,k)  +ua(j,i,k))*twt(k,1)) *       &
                  (f(j,i,k)+f(j-1,i,k)) +                     &
                 ((va(j,i+1,k-1)+va(j+1,i+1,k-1))*twt(k,2) +  &
                  (va(j,i+1,k)  +va(j+1,i+1,k))*twt(k,1)) *   &
                  (f(j,i+1,k)+f(j,i,k)) -                     &
                 ((va(j,i,k-1)+va(j+1,i,k-1))*twt(k,2) +      &
                  (va(j,i,k)  +va(j+1,i,k))*twt(k,1)) *       &
                  (f(j,i-1,k)+f(j,i,k)))
            end do
          end do
        end do
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine hadv3d

    subroutine hadv4d(ften,f,nk,m,p)
      implicit none
      integer(ik4) , intent (in) :: nk , m , p
      real(rk8) , pointer , intent (in) , dimension(:,:,:,:) :: f
      real(rk8) , pointer , intent (inout), dimension(:,:,:,:) :: ften

      integer(ik4) :: i , j , k , n
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'hadv4d'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      !
      ! for qv:
      !
      do n = m , p
        do k = 1 , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,k,n) = ften(j,i,k,n) - xmapf(j,i)       *            &
                  ((ua(j+1,i+1,k)+ua(j+1,i,k))*(f(j+1,i,k,n)+f(j,i,k,n)) -  &
                   (ua(j,i+1,k)+ua(j,i,k)) *   (f(j,i,k,n)+f(j-1,i,k,n)) +  &
                   (va(j+1,i+1,k)+va(j,i+1,k))*(f(j,i+1,k,n)+f(j,i,k,n)) -  &
                   (va(j+1,i,k)+va(j,i,k)) *   (f(j,i-1,k,n)+f(j,i,k,n)))
            end do
          end do
        end do
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine hadv4d
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
      integer(ik4) , intent(in) :: ind , nk
      real(rk8) , pointer , intent (in) , dimension(:,:,:) :: f
      real(rk8) , pointer , intent (inout), dimension(:,:,:) :: ften

      real(rk8) :: slope , rdphf , rdplf , ff , qq
      integer(ik4) :: i , j , k
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'vadv3d'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      if ( ind == 0 ) then
        if ( nk == kz ) then
          do k = 2 , nk
            do i = ici1 , ici2
              do j = jci1 , jci2
                ff = (twt(k,1)*f(j,i,k)+twt(k,2)*f(j,i,k-1)) * svv(j,i,k)
                ften(j,i,k-1) = ften(j,i,k-1) - ff*xds(k-1)
                ften(j,i,k)   = ften(j,i,k)   + ff*xds(k)
              end do
            end do
          end do
        else
          do k = 1 , nk - 1
            do i = ici1 , ici2
              do j = jci1 , jci2
                qq = d_half * (svv(j,i,k) + svv(j,i,k+1))
                ff = qq * (f(j,i,k) + f(j,i,k+1))
                ften(j,i,k+1) = ften(j,i,k+1) + ff*dds(k+1)
                ften(j,i,k)   = ften(j,i,k)   - ff*dds(k)
              end do
            end do
          end do
        end if
      else if ( ind == 1 ) then
        !
        ! vertical advection terms : interpolate to full sigma levels
        !
        do k = 2 , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              rdphf = (pfs(j,i,k)/phs(j,i,k))**c287
              rdplf = (pfs(j,i,k)/phs(j,i,k-1))**c287
              ff = (twt(k,1)*f(j,i,k)  *rdphf + &
                    twt(k,2)*f(j,i,k-1)*rdplf) * svv(j,i,k)
              ften(j,i,k-1) = ften(j,i,k-1) - ff*xds(k-1)
              ften(j,i,k)   = ften(j,i,k)   + ff*xds(k)
            end do
          end do
        end do
      else if ( ind == 2 ) then
        fg(:,:,1) = d_zero
        do k = 2 , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              fg(j,i,k)= twt(k,1)*f(j,i,k) + twt(k,2)*f(j,i,k-1)
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( kpb(j,i) > nk ) then
              call fatal(__FILE__,__LINE__,'kpbl is greater than nk')
            end if
            if ( kpb(j,i) >= 4 ) then
              ! Calculate slope of scalar in layer above ambiguous layer
              k = kpb(j,i)-2
              if ( (f(j,i,k+1)-f(j,i,k)) > d_zero .and. &
                   (f(j,i,k)-f(j,i,k-1)) > d_zero ) then
                slope = min((f(j,i,k+1)-f(j,i,k))/(hsigma(k+1)-hsigma(k)), &
                            (f(j,i,k)-f(j,i,k-1))/(hsigma(k)-hsigma(k-1)))
              else if ( (f(j,i,k+1)-f(j,i,k)) < d_zero .and. &
                        (f(j,i,k)-f(j,i,k-1)) < d_zero ) then
                slope = max((f(j,i,k+1)-f(j,i,k))/(hsigma(k+1)-hsigma(k)), &
                            (f(j,i,k)-f(j,i,k-1))/(hsigma(k)-hsigma(k-1)))
              else
                slope = d_zero
              end if
              ! Replace the values of scalar at top and bottom of ambiguous
              ! layer as long as inversion is actually in the ambiguous layer
              k = kpb(j,i)
              fg(j,i,k-1) = f(j,i,k-2) + slope*(sigma(k-1)-hsigma(k-2))
              if (abs(f(j,i,k-2) + slope*(hsigma(k-1)-hsigma(k-2))-f(j,i,k)) > &
                  abs(f(j,i,k-1)-f(j,i,k)) ) then
                fg(j,i,k) = f(j,i,k)
              else
                fg(j,i,k) = f(j,i,k-2) + slope*(sigma(k)-hsigma(k-2))
              end if
            end if
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            ften(j,i,1) = ften(j,i,1) - svv(j,i,2)*fg(j,i,2)*xds(1)
          end do
        end do
        do k = 2 , nk-1
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,k) = ften(j,i,k) - &
                  (svv(j,i,k+1)*fg(j,i,k+1)-svv(j,i,k)*fg(j,i,k))*xds(k)
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            ften(j,i,nk) = ften(j,i,nk) + &
                   svv(j,i,nk)*fg(j,i,nk)*xds(nk)
          end do
        end do
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine vadv3d
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
    subroutine vadv4d(ften,f,nk,ind,m,p)
      implicit none
      integer(ik4) , intent(in) :: ind , nk
      integer(ik4) , optional , intent(in) :: m , p
      real(rk8) , pointer , intent (in) , dimension(:,:,:,:) :: f
      real(rk8) , pointer , intent (inout), dimension(:,:,:,:) :: ften

      real(rk8) :: slope , ff
      integer(ik4) :: i , j , k , n , n1 , n2
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'vadv4d'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      if ( present(m) ) then
        if ( present(p) ) then
          n1 = m
          n2 = p
        else
          n1 = m
          n2 = m
        end if
      else
        n1 = lbound(f,4)
        n2 = ubound(f,4)
      end if
      if ( ind == 1 ) then
        do n = n1 , n2
          do k = 2 , nk
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( f(j,i,k,n) > minqv .and. f(j,i,k-1,n) > minqv ) then
                  ff = (f(j,i,k,n) * &
                    (f(j,i,k-1,n)/f(j,i,k,n))**qcon(k)) * svv(j,i,k)
                else
                  ff = d_zero
                end if
                ften(j,i,k-1,n) = ften(j,i,k-1,n) - ff*xds(k-1)
                ften(j,i,k,n)   = ften(j,i,k,n)   + ff*xds(k)
              end do
            end do
          end do
        end do
      else if ( ind == 2 ) then
        do n = n1 , n2
          do k = 2 , nk
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( f(j,i,k,n) > mintr .and. f(j,i,k-1,n) > mintr ) then
                  ff = (twt(k,1)*f(j,i,k,n)+twt(k,2)*f(j,i,k-1,n)) * svv(j,i,k)
                else
                  ff = d_zero
                end if
                ften(j,i,k-1,n) = ften(j,i,k-1,n) - ff*xds(k-1)
                ften(j,i,k,n)   = ften(j,i,k,n)   + ff*xds(k)
              end do
            end do
          end do
        end do
      else if ( ind == 3 ) then
        do n = n1 , n2
          do k = 2 , nk
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( f(j,i,k,n) > minqx .and. f(j,i,k-1,n) > minqx ) then
                  ff = (twt(k,1)*f(j,i,k,n)+twt(k,2)*f(j,i,k-1,n)) * svv(j,i,k)
                else
                  ff = d_zero
                end if
                ften(j,i,k-1,n) = ften(j,i,k-1,n) - ff*xds(k-1)
                ften(j,i,k,n)   = ften(j,i,k,n)   + ff*xds(k)
              end do
            end do
          end do
        end do
      else if ( ind == 4 ) then
        do n = n1 , n2
          fg(:,:,1) = d_zero
          do k = 2 , nk
            do i = ici1 , ici2
              do j = jci1 , jci2
                fg(j,i,k) = twt(k,1)*f(j,i,k,n) + twt(k,2)*f(j,i,k-1,n)
              end do
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( kpb(j,i) > nk ) then
                call fatal(__FILE__,__LINE__,'kpbl is greater than nk')
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
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,1,n) = ften(j,i,1,n) - svv(j,i,2)*fg(j,i,2)*xds(1)
            end do
          end do
          do k = 2 , nk-1
            do i = ici1 , ici2
              do j = jci1 , jci2
                ften(j,i,k,n) = ften(j,i,k,n) - &
                    (svv(j,i,k+1)*fg(j,i,k+1)-svv(j,i,k)*fg(j,i,k))*xds(k)
              end do
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,nk,n) = ften(j,i,nk,n) + &
                     svv(j,i,nk)*fg(j,i,nk)*xds(nk)
            end do
          end do
        end do
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine vadv4d

end module mod_advection

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
