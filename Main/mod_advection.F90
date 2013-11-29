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
  use mod_service

  private
   
  public :: init_advection, hadv , vadv

  interface hadv
    module procedure hadv3d
    module procedure hadv4d
  end interface hadv

  interface vadv
    module procedure vadv3d
    module procedure vadv4d
  end interface vadv

  real(rk8) , pointer , dimension(:,:,:) :: ua   ! U wind * ps
  real(rk8) , pointer , dimension(:,:,:) :: va   ! V wind * ps
  real(rk8) , pointer , dimension(:,:) :: ps     ! Surface pressure
  real(rk8) , pointer , dimension(:,:) :: mapfx  ! Map factor Cross
  real(rk8) , pointer , dimension(:,:) :: mapfd  ! Map factor Dot
  real(rk8) , pointer , dimension(:,:,:) :: vsv  ! Vertical Sigma Velocity
  integer(ik4) , pointer , dimension(:,:) :: kpbl   ! Top of PBL
!
! working space used to store the interlated values in vadv.
!
  real(rk8) , pointer , dimension(:,:,:) :: fg

  real(rk8) , parameter :: c287 = 0.287D+00
!
  real(rk8) , parameter :: falow = 1.0D-8
!
  contains

    subroutine init_advection(dom,sfs,atm,vertvel,kpbltop)
      implicit none
      type(domain) , intent(in) :: dom
      type(surfstate), intent(in) :: sfs
      type(atmstate) , intent(in) :: atm
      real(rk8) , pointer , dimension(:,:,:) :: vertvel
      integer(ik4) , pointer , dimension(:,:) :: kpbltop
      call assignpnt(atm%u,ua)
      call assignpnt(atm%v,va)
      call assignpnt(sfs%psa,ps)
      call assignpnt(dom%msfx,mapfx)
      call assignpnt(dom%msfd,mapfd)
      call assignpnt(vertvel,vsv)
      call assignpnt(kpbltop,kpbl)
      call getmem3d(fg,jde1,jde2,ide1,ide2,1,kz,'mod_advection:fg')
    end subroutine init_advection
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!  HADV                                                               c
!                                                                     c
!     This subroutines computes the horizontal flux-divergence terms. c
!     second-order difference is used.                                c
!                                                                     c
!     ldot   : cross/dot variable flag                                c
!                                                                     c
!     ften   : is the tendency for variable 'f'.                      c
!                                                                     c
!     f      : is p*f.                                                c
!                                                                     c
!     nk     : is the number of vertical levels to work (kz/kzp1)     c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine hadv3d(ldot,ften,f,nk)
      implicit none
      logical , intent(in) :: ldot ! Cross/dot flag
      integer(ik4) , intent (in) :: nk
      real(rk8) , pointer , intent (in) , dimension(:,:,:) :: f
      real(rk8) , pointer , intent (inout), dimension(:,:,:) :: ften
!
      real(rk8) :: ucmona , ucmonb , ucmonc , vcmona , vcmonb , vcmonc
      integer(ik4) :: i , j , k
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'hadv3d'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      if ( ldot ) then
        !
        ! ua, va : are p*u and p*v.
        ! msfd   : is the map scale factor at dot points.
        !
        do k = 1 , nk
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
              ften(j,i,k) = ften(j,i,k) -                  &
                          ((f(j+1,i,k)+f(j,i,k))*ucmonb -  &
                           (f(j,i,k)+f(j-1,i,k))*ucmonc +  &
                           (f(j,i+1,k)+f(j,i,k))*vcmonb -  &
                           (f(j,i,k)+f(j,i-1,k))*vcmonc) / &
                           (dx16*mapfd(j,i)*mapfd(j,i))
            end do
          end do
        end do
      else
        !
        ! for t
        !
        do k = 1 , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,k) = ften(j,i,k) -                               &
                  ((ua(j+1,i+1,k)+ua(j+1,i,k))*(f(j+1,i,k)+f(j,i,k)) -  &
                   (ua(j,i+1,k)+ua(j,i,k)) *   (f(j,i,k)+f(j-1,i,k)) +  &
                   (va(j+1,i+1,k)+va(j,i+1,k))*(f(j,i+1,k)+f(j,i,k)) -  &
                   (va(j+1,i,k)+va(j,i,k)) *   (f(j,i-1,k)+f(j,i,k))) / &
                   (dx4*mapfx(j,i)*mapfx(j,i))
            end do
          end do
        end do
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine hadv3d
!
    subroutine hadv4d(ften,f,nk,m,p)
      implicit none
      integer(ik4) , intent (in) :: nk
      integer(ik4) , optional , intent (in) :: m , p
      real(rk8) , pointer , intent (in) , dimension(:,:,:,:) :: f
      real(rk8) , pointer , intent (inout), dimension(:,:,:,:) :: ften
!
      integer(ik4) :: i , j , k , n , n1 , n2
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'hadv4d'
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
      !
      ! for qv:
      !
      do n = n1 , n2
        do k = 1 , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,k,n) = ften(j,i,k,n) -                               &
                  ((ua(j+1,i+1,k)+ua(j+1,i,k))*(f(j+1,i,k,n)+f(j,i,k,n)) -  &
                   (ua(j,i+1,k)+ua(j,i,k)) *   (f(j,i,k,n)+f(j-1,i,k,n)) +  &
                   (va(j+1,i+1,k)+va(j,i+1,k))*(f(j,i+1,k,n)+f(j,i,k,n)) -  &
                   (va(j+1,i,k)+va(j,i,k)) *   (f(j,i-1,k,n)+f(j,i,k,n))) / &
                   (dx4*mapfx(j,i)*mapfx(j,i))
            end do
          end do
        end do
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine hadv4d
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
! VADV                                                                c
!                                                                     c
!     This subroutine computes the vertical flux-divergence terms.    c
!                                                                     c
!     ften   : is the tendency of variable 'f'.                       c
!                                                                     c
!     f      : is p*f.                                                c
!                                                                     c
!     jstart : is the j'th slice of f anf ften to start               c
!                                                                     c
!     jsstop : is the j'th slice of f anf ften to stop                c
!                                                                     c
!     ind = 1 : for t.                                                c
!           2 : for qv.                                               c
!           3 : for qc and qr.                                        c
!           4 : for u and v.                                          c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine vadv3d(ldot,ften,f,nk,ind)
      implicit none
      logical , intent(in) :: ldot
      integer(ik4) , intent(in) :: ind , nk
      real(rk8) , pointer , intent (in) , dimension(:,:,:) :: f
      real(rk8) , pointer , intent (inout), dimension(:,:,:) :: ften
!
      real(rk8) :: slope
      integer(ik4) :: i , j , k
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'vadv3d'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      if ( ldot ) then
        if ( ind /= 2 ) then
          call fatal(__FILE__,__LINE__, &
                     'The advection scheme you required is not available.')
        end if
      end if

      if ( ind == 1 ) then
        !
        ! vertical advection terms : interpolate t to full sigma levels
        !
        do i = ici1 , ici2
          do j = jci1 , jci2
            fg(j,i,1) = d_zero
          end do
        end do
        do k = 2 , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              fg(j,i,k) = twt(k,1)*f(j,i,k) *                                 &
                   ((ps(j,i)*sigma(k)+ptop)/(ps(j,i)*hsigma(k)+ptop))**c287 + &
                        twt(k,2)*f(j,i,k-1) *                                 &
                   ((ps(j,i)*sigma(k)+ptop)/(ps(j,i)*hsigma(k-1)+ptop))**c287
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            ften(j,i,1) = ften(j,i,1) - vsv(j,i,2)*fg(j,i,2)/dsigma(1)
          end do
        end do
        do k = 2 , nk-1
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,k) = ften(j,i,k) - &
                   (vsv(j,i,k+1)*fg(j,i,k+1)-vsv(j,i,k)*fg(j,i,k))/dsigma(k)
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            ften(j,i,nk) = ften(j,i,nk) + vsv(j,i,nk)*fg(j,i,nk)/dsigma(nk)
          end do
        end do
      else if ( ind == 2 ) then
        !
        ! vertical advection terms : interpolate ua or va to full sigma levels
        !
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            fg(j,i,1) = d_zero
          end do
        end do
        do k = 2 , nk
          do i = idi1 , idi2
            do j = jdi1 , jdi2
              fg(j,i,k) = d_half*(f(j,i,k)+f(j,i,k-1))/mapfd(j,i)
            end do
          end do
        end do
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            ften(j,i,1) = ften(j,i,1) -                &
                        (vsv(j-1,i-1,2)+vsv(j-1,i,2) + &
                         vsv(j,i,2)+vsv(j,i-1,2))    * &
                         fg(j,i,2)/(d_four*dsigma(1))
          end do
        end do
        do k = 2 , nk-1
          do i = idi1 , idi2
            do j = jdi1 , jdi2
              ften(j,i,k) = ften(j,i,k) -                                &
                          ((vsv(j-1,i,k+1)+vsv(j-1,i-1,k+1)+             &
                            vsv(j,i,k+1)  +vsv(j,i-1,k+1))*fg(j,i,k+1) - &
                           (vsv(j-1,i,k)  +vsv(j-1,i-1,k)+               &
                            vsv(j,i,k)    +vsv(j,i-1,k))*fg(j,i,k)) /    &
                          (d_four*dsigma(k))
            end do
          end do
        end do
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            ften(j,i,nk) = ften(j,i,nk) +                 &
                         (vsv(j-1,i,nk)+vsv(j-1,i-1,nk) + &
                          vsv(j,i,nk)+vsv(j,i-1,nk)) *    &
                          fg(j,i,nk)/(d_four*dsigma(nk))
          end do
        end do
      else if ( ind == 3 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            fg(j,i,1) = d_zero
          end do
        end do
        do k = 2 , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              fg(j,i,k)= twt(k,1)*f(j,i,k) + twt(k,2)*f(j,i,k-1)
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( kpbl(j,i).gt.nk ) then
              call fatal(__FILE__,__LINE__,'kpbl is greater than nk')
            end if
            if ( kpbl(j,i).ge.4 ) then
              ! Calculate slope of scalar in layer above ambiguous layer
              k = kpbl(j,i)-2
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
              k = kpbl(j,i)
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
            ften(j,i,1) = ften(j,i,1) - vsv(j,i,2)*fg(j,i,2)/dsigma(1)
          end do
        end do
        do k = 2 , nk-1
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,k) = ften(j,i,k) - &
                  (vsv(j,i,k+1)*fg(j,i,k+1)-vsv(j,i,k)*fg(j,i,k))/dsigma(k)
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            ften(j,i,nk) = ften(j,i,nk) + &
                   vsv(j,i,nk)*fg(j,i,nk)/dsigma(nk)
          end do
        end do
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine vadv3d
!
    subroutine vadv4d(ften,f,nk,ind,m,p)
      implicit none
      integer(ik4) , intent(in) :: ind , nk
      integer(ik4) , optional , intent(in) :: m , p
      real(rk8) , pointer , intent (in) , dimension(:,:,:,:) :: f
      real(rk8) , pointer , intent (inout), dimension(:,:,:,:) :: ften
!
      real(rk8) :: slope
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
          !
          ! vertical advection term : interpolate qv to full sigma levels
          !
          do i = ici1 , ici2
            do j = jci1 , jci2
              fg(j,i,1) = d_zero
            end do
          end do
          do k = 2 , nk
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( f(j,i,k,n) > falow .and. f(j,i,k-1,n) > falow ) then
                  fg(j,i,k) = f(j,i,k,n)*(f(j,i,k-1,n)/f(j,i,k,n))**qcon(k)
                else
                  fg(j,i,k) = d_zero
                end if
              end do
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,1,n) = ften(j,i,1,n) - vsv(j,i,2)*fg(j,i,2)/dsigma(1)
            end do
          end do
          do k = 2 , nk-1
            do i = ici1 , ici2
              do j = jci1 , jci2
                ften(j,i,k,n) = ften(j,i,k,n) - &
                       (vsv(j,i,k+1)*fg(j,i,k+1)-vsv(j,i,k)*fg(j,i,k))/dsigma(k)
              end do
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,nk,n) = ften(j,i,nk,n)+vsv(j,i,nk)*fg(j,i,nk)/dsigma(nk)
            end do
          end do
        end do
      else if ( ind == 2 ) then
        do n = n1 , n2
          do i = ici1 , ici2
            do j = jci1 , jci2
              fg(j,i,1) = d_zero
            end do
          end do
          do k = 2 , nk
            do i = ici1 , ici2
              do j = jci1 , jci2
                fg(j,i,k) = twt(k,1)*f(j,i,k,n) + twt(k,2)*f(j,i,k-1,n)
              end do
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,1,n) = ften(j,i,1,n) - vsv(j,i,2)*fg(j,i,2)/dsigma(1)
            end do
          end do
          do k = 2 , nk-1
            do i = ici1 , ici2
              do j = jci1 , jci2
                ften(j,i,k,n) = ften(j,i,k,n) - &
                       (vsv(j,i,k+1)*fg(j,i,k+1)-vsv(j,i,k)*fg(j,i,k))/dsigma(k)
              end do
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,nk,n) = ften(j,i,nk,n)+vsv(j,i,nk)*fg(j,i,nk)/dsigma(nk)
            end do
          end do
        end do
      else if ( ind == 3 ) then
        do n = n1 , n2
          do i = ici1 , ici2
            do j = jci1 , jci2
              fg(j,i,1) = d_zero
            end do
          end do
          do k = 2 , nk
            do i = ici1 , ici2
              do j = jci1 , jci2
                fg(j,i,k)= twt(k,1)*f(j,i,k,n) + twt(k,2)*f(j,i,k-1,n)
              end do
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( kpbl(j,i).gt.nk ) then
                call fatal(__FILE__,__LINE__,'kpbl is greater than nk')
              end if
              if ( kpbl(j,i).ge.4 ) then
                ! Calculate slope of scalar in layer above ambiguous layer
                k = kpbl(j,i)-2
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
                k = kpbl(j,i)
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
              ften(j,i,1,n) = ften(j,i,1,n) - vsv(j,i,2)*fg(j,i,2)/dsigma(1)
            end do
          end do
          do k = 2 , nk-1
            do i = ici1 , ici2
              do j = jci1 , jci2
                ften(j,i,k,n) = ften(j,i,k,n) - &
                    (vsv(j,i,k+1)*fg(j,i,k+1)-vsv(j,i,k)*fg(j,i,k))/dsigma(k)
              end do
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              ften(j,i,nk,n) = ften(j,i,nk,n) + &
                     vsv(j,i,nk)*fg(j,i,nk)/dsigma(nk)
            end do
          end do
        end do
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine vadv4d
!
end module mod_advection
