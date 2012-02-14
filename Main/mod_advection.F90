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
  use mod_atm_interface
  use mod_dynparam
  use mod_runparams
  use mod_memutil
  use mod_service

  private
   
  public :: init_advection, hadv , vadv

  real(dp) , pointer , dimension(:,:,:) :: u    ! U wind * ps
  real(dp) , pointer , dimension(:,:,:) :: v    ! V wind * ps
  real(dp) , pointer , dimension(:,:) :: ps     ! Surface pressure
  real(dp) , pointer , dimension(:,:) :: mapfx  ! Map factor Cross
  real(dp) , pointer , dimension(:,:) :: mapfd  ! Map factor Dot
  real(dp) , pointer , dimension(:,:,:) :: vsv  ! Vertical Sigma Velocity
  integer , pointer , dimension(:,:) :: kpbl   ! Top of PBL
!
! working space used to store the interlated values in vadv.
!
  real(dp) , pointer , dimension(:,:,:) :: fg

  real(dp) , parameter :: c287 = 0.287D+00
!
! relaxed upstream scheme factors
!
  real(dp) , parameter :: fact1 = 0.60D0
!hy
! real(dp) , parameter :: fact1 = 0.75D0
!hy
  real(dp) , parameter :: fact2 = d_one - fact1
  real(dp) , parameter :: falow = 1.0D-8
!
  contains

    subroutine init_advection(dom,sfs,atm,vertvel,kpbltop)
      implicit none
      type(domain) , intent(in) :: dom
      type(surfstate), intent(in) :: sfs
      type(atmstate) , intent(in) :: atm
      real(dp) , pointer , dimension(:,:,:) :: vertvel
      integer , pointer , dimension(:,:) :: kpbltop

      call assignpnt(atm%u,u)
      call assignpnt(atm%v,v)
      call assignpnt(sfs%psa,ps)
      call assignpnt(dom%msfx,mapfx)
      call assignpnt(dom%msfd,mapfd)
      call assignpnt(vertvel,vsv)
      call assignpnt(kpbltop,kpbl)
      call getmem3d(fg,1,jxp,idot1,idot2,1,kz,'mod_advection:fg')
    end subroutine init_advection
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!  HADV                                                               c
!                                                                     c
!     This subroutines computes the horizontal flux-divergence terms. c
!     second-order difference is used.                                c
!                                                                     c
!     ldot   : cross/dot variable flagg                               c
!                                                                     c
!     ften   : is the tendency for variable 'f'.                      c
!                                                                     c
!     f      : is p*f.                                                c
!                                                                     c
!     jstart : is the j'th slice of f anf ften to start               c
!                                                                     c
!     jsstop : is the j'th slice of f anf ften to stop                c
!                                                                     c
!     ind = 1 : for t and qv.                                         c
!         = 2 : for qc and qr.                                        c
!         = 3 : for u and v.                                          c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine hadv(ldot,ften,f,nk,ind)
!
      implicit none
!
      logical , intent(in) :: ldot ! Cross/dot flag
      integer , intent (in) :: ind , nk
      real(dp) , pointer , intent (in) , dimension(:,:,:) :: f
      real(dp) , pointer , intent (inout), dimension(:,:,:) :: ften
!
      real(dp) :: fx1 , fx2 , fy1 , fy2
      real(dp) :: ucmona , ucmonb , ucmonc , vcmona , vcmonb , vcmonc
      integer :: i , j , k
!
      character (len=64) :: subroutine_name='hadv'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
      if ( ldot ) then
        if ( ind == 3 ) then
          !
          ! ua, va : are p*u and p*v.
          ! msfd   : is the map scale factor at dot points.
          !
          do k = 1 , nk
            do i = idi1 , idi2
              do j = jdi1 , jdi2
                ucmona = u(j,i+1,k)+d_two*u(j,i,k)+u(j,i-1,k)
                vcmona = v(j+1,i,k)+d_two*v(j,i,k)+v(j-1,i,k)
                ucmonb = u(j+1,i+1,k) + d_two*u(j+1,i,k) + &
                         u(j+1,i-1,k) + ucmona
                vcmonb = v(j+1,i+1,k) + d_two*v(j,i+1,k) + &
                         v(j-1,i+1,k) + vcmona
                ucmonc = u(j-1,i+1,k) + d_two*u(j-1,i,k) + &
                         u(j-1,i-1,k) + ucmona
                vcmonc = v(j+1,i-1,k) + d_two*v(j,i-1,k) + &
                         v(j-1,i-1,k) + vcmona
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
          call fatal(__FILE__,__LINE__, &
                     'The advection scheme you required is not available.')
        end if
      else
        if ( ind == 1 ) then
          !
          ! for t and qv:
          !
          do k = 1 , nk
            do i = ici1 , ici2
              do j = jci1 , jci2
                ften(j,i,k) = ften(j,i,k) -                             &
                    ((u(j+1,i+1,k)+u(j+1,i,k))*(f(j+1,i,k)+f(j,i,k)) -  &
                     (u(j,i+1,k)+u(j,i,k)) *   (f(j,i,k)+f(j-1,i,k)) +  &
                     (v(j+1,i+1,k)+v(j,i+1,k))*(f(j,i+1,k)+f(j,i,k)) -  &
                     (v(j+1,i,k)+v(j,i,k)) *   (f(j,i-1,k)+f(j,i,k))) / &
                     (dx4*mapfx(j,i)*mapfx(j,i))
              end do
            end do
          end do
        else if ( ind == 2 ) then
          !
          ! for qc and qr:
          ! up-wind values of qc and qr are used.
          !
          do k = 1 , nk
            do i = ici1 , ici2
              do j = jci1 , jci2
                ucmonb = d_half*(u(j+1,i+1,k)+u(j+1,i,k))
                ucmona = d_half*(u(j,i+1,k)+u(j,i,k))
                if ( ucmonb >= d_zero ) then
                  fx2 = fact1*f(j,i,k) + fact2*f(j+1,i,k)
                else
                  fx2 = fact1*f(j+1,i,k) + fact2*f(j,i,k)
                end if
                if ( ucmona >= d_zero ) then
                  fx1 = fact1*f(j-1,i,k) + fact2*f(j,i,k)
                else
                  fx1 = fact1*f(j,i,k) + fact2*f(j-1,i,k)
                end if
                vcmonb = d_half*(v(j+1,i+1,k)+v(j,i+1,k))
                vcmona = d_half*(v(j+1,i,k)+v(j,i,k))
                if ( vcmonb >= d_zero ) then
                  fy2 = fact1*f(j,i,k) + fact2*f(j,i+1,k)
                else
                  fy2 = fact1*f(j,i+1,k) + fact2*f(j,i,k)
                end if
                if ( vcmona >= d_zero ) then
                  fy1 = fact1*f(j,i+1,k) + fact2*f(j,i,k)
                else
                  fy1 = fact1*f(j,i,k) + fact2*f(j,i-1,k)
                end if
                ften(j,i,k) = ften(j,i,k) -                          &
                     (ucmonb*fx2-ucmona*fx1+vcmonb*fy2-vcmona*fy1) / &
                     (dx*mapfx(j,i)*mapfx(j,i))
              end do
            end do
          end do
        else
          call fatal(__FILE__,__LINE__, &
                     'The advection scheme you required is not available.')
        end if
      end if
      call time_end(subroutine_name,idindx)
!
    end subroutine hadv
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
    subroutine vadv(ldot,ften,f,nk,ind)
!
      implicit none
!
      logical , intent(in) :: ldot
      integer , intent(in) :: ind , nk
      real(dp) , pointer , intent (in) , dimension(:,:,:) :: f
      real(dp) , pointer , intent (inout), dimension(:,:,:) :: ften
!
      real(dp) :: f1 , f2 , slope
      integer :: i , j , k
!
      character (len=64) :: subroutine_name='vadv'
      integer :: idindx=0
!
!----------------------------------------------------------------------
!
      call time_begin(subroutine_name,idindx)
!
      if ( ldot ) then
        if ( ind /= 4 ) then
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
              fg(j,i,k) = twt(k,1)*f(j,i,k) *                              &
                   ((ps(j,i)*sigma(k)+ptop)/(ps(j,i)*a(k)+ptop))**c287 + &
                        twt(k,2)*f(j,i,k-1) *                            &
                   ((ps(j,i)*sigma(k)+ptop)/(ps(j,i)*a(k-1)+ptop))**c287
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
!
      else if ( ind == 2 ) then
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
              if ( f(j,i,k) > falow .and. f(j,i,k-1) > falow ) then
                fg(j,i,k) = f(j,i,k)*(f(j,i,k-1)/f(j,i,k))**qcon(k)
              else
                fg(j,i,k) = d_zero
              end if
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            ften(j,i,1) = ften(j,i,1) - &
                  vsv(j,i,2)*fg(j,i,2)/dsigma(1)
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
!
      else if ( ind == 3 ) then
        !
        ! vertical advection terms for qc and qr:
        !
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( vsv(j,i,2) >= d_zero ) then
              f2 = f(j,i,1)
            else
              f2 = f(j,i,2)
            end if
            ften(j,i,1) = ften(j,i,1) - vsv(j,i,2)*f2/dsigma(1)
          end do
        end do
        do k = 2 , nk-1
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( vsv(j,i,k+1) >= d_zero ) then
                f2 = f(j,i,k)
              else
                f2 = f(j,i,k+1)
              end if
              if ( vsv(j,i,k) >= d_zero ) then
                f1 = f(j,i,k-1)
              else
                f1 = f(j,i,k)
              end if
              ften(j,i,k) = ften(j,i,k) - &
                    (vsv(j,i,k+1)*f2-vsv(j,i,k)*f1)/dsigma(k)
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( vsv(j,i,nk) >= d_zero ) then
              f1 = f(j,i,nk-1)
            else
              f1 = f(j,i,nk)
            end if
            ften(j,i,nk) = ften(j,i,nk) + vsv(j,i,nk)*f1/dsigma(nk)
          end do
        end do
!
      else if ( ind == 4 ) then
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
!
      else if ( ind == 5 ) then
        do k = 2 , nk
          do i = ici1 , ici2
            do j = jci1 , jci2
              fg(j,i,k) = twt(k,1)*f(j,i,k) + twt(k,2)*f(j,i,k-1)
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
!
      else if ( ind == 6 ) then
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
                slope = min((f(j,i,k+1)-f(j,i,k))/(a(k+1)-a(k)), &
                            (f(j,i,k)-f(j,i,k-1))/(a(k)-a(k-1)))
              else if ( (f(j,i,k+1)-f(j,i,k)) < d_zero .and. &
                        (f(j,i,k)-f(j,i,k-1)) < d_zero ) then
                slope = max((f(j,i,k+1)-f(j,i,k))/(a(k+1)-a(k)), &
                            (f(j,i,k)-f(j,i,k-1))/(a(k)-a(k-1)))
              else
                slope = d_zero
              end if
              ! Replace the values of scalar at top and bottom of ambiguous
              ! layer as long as inversion is actually in the ambiguous layer
              k = kpbl(j,i)
              fg(j,i,k-1) = f(j,i,k-2) + slope*(sigma(k-1)-a(k-2))
              if (abs(f(j,i,k-2) + slope*(a(k-1)-a(k-2))-f(j,i,k)) > &
                  abs(f(j,i,k-1)-f(j,i,k)) ) then
                fg(j,i,k) = f(j,i,k)
              else
                fg(j,i,k) = f(j,i,k-2) + slope*(sigma(k)-a(k-2))
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
!
      call time_end(subroutine_name,idindx)

    end subroutine vadv
!
end module mod_advection
