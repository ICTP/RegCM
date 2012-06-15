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
 
module mod_diffusion
!
! Diffusion calculations
!
  use mod_atm_interface
  use mod_runparams
  use mod_mppparam
  use mod_service 
  private

  interface diffu_x
    module procedure diffu_x3d
    module procedure diffu_x4d
  end interface

  public :: diffu_d , diffu_x
!
  contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     These subroutines computes the diffusion term for decoupled     c
!     variable on constant sigma surface.                             c
!                                                                     c
!     ften    : tendency for variable                                 c
!                                                                     c
!     xkc     : horizontal diffusion coefficient                      c
!                                                                     c
!     f       : variable f on cross/dot points                        c
!                                                                     c
!     ind = 1 : var is already multiplied by map scale factor         c
!         = 0 : var is "not"   multiplied by map scale factor         c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine diffu_d(ften,f,press,mapf,xkc,ind)
!
    implicit none
!
    integer , intent(in) :: ind
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: f
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: xkc
    real(dp) , pointer , dimension(:,:) , intent(in) :: press
    real(dp) , pointer , dimension(:,:) , intent(in) :: mapf
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ften
!
    integer :: i , j , k
!
    character (len=64) :: subroutine_name='diffu_d'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)

    !
    ! fourth-order scheme for interior:
    !
    do k = 1 , kz
      do i = idii1 , idii2
        do j = jdii1 , jdii2
          if ( ind == 0 ) then
            ften(j,i,k) = ften(j,i,k) - xkc(j,i,k) *      &
                  rdxsq*(f(j+2,i,k)+f(j-2,i,k)+f(j,i+2,k)+f(j,i-2,k)  - &
                 d_four*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) + &
                d_twelve*f(j,i,k))*press(j,i)
          else
            ften(j,i,k) = ften(j,i,k) - xkc(j,i,k) *        &
                  rdxsq*(f(j+2,i,k)/mapf(j+2,i)+f(j-2,i,k)/mapf(j-2,i) +  &
                         f(j,i+2,k)/mapf(j,i+2)+f(j,i-2,k)/mapf(j,i-2) -  &
                 d_four*(f(j+1,i,k)/mapf(j+1,i)+f(j-1,i,k)/mapf(j-1,i) +  &
                         f(j,i+1,k)/mapf(j,i+1)+f(j,i-1,k)/mapf(j,i-1)) + &
                d_twelve*f(j,i,k)/mapf(j,i))*press(j,i)
          end if
        end do
      end do
    end do
    !
    ! second-order scheme for east and west boundaries:
    !
    if ( ma%hasleft ) then
      j = jdi1
      do k = 1 , kz
        do i = idi1 , idi2
          if ( ind == 0 ) then
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *    &
                          rdxsq*(f(j+1,i,k)+f(j-1,i,k) + &
                                 f(j,i+1,k)+f(j,i-1,k) - &
                          d_four*f(j,i,k))*press(j,i)
          else
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *     &
                          rdxsq*(f(j+1,i,k)/mapf(j+1,i) + &
                                 f(j-1,i,k)/mapf(j-1,i) + &
                                 f(j,i+1,k)/mapf(j,i+1) + &
                                 f(j,i-1,k)/mapf(j,i-1) - &
                          d_four*f(j,i,k)/mapf(j,i))*press(j,i)
          end if
        end do
      end do
    end if
    if ( ma%hasright ) then
      j = jdi2
      do k = 1 , kz
        do i = idi1 , idi2
          if ( ind == 0 ) then
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *    &
                          rdxsq*(f(j+1,i,k)+f(j-1,i,k) + &
                                 f(j,i+1,k)+f(j,i-1,k) - &
                          d_four*f(j,i,k))*press(j,i)
          else
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *     &
                          rdxsq*(f(j+1,i,k)/mapf(j+1,i) + &
                                 f(j-1,i,k)/mapf(j-1,i) + &
                                 f(j,i+1,k)/mapf(j,i+1) + &
                                 f(j,i-1,k)/mapf(j,i-1) - &
                          d_four*f(j,i,k)/mapf(j,i))*press(j,i)
          end if
        end do
      end do
    end if
    !
    ! second-order scheme for north and south boundaries:
    !
    if ( ma%hasbottom ) then
      i = idi1
      do k = 1 , kz
        do j = jdi1 , jdi2
          if ( ind == 0 ) then
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *         &
                  rdxsq*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k) - &
                  d_four*f(j,i,k))*press(j,i)
          else
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *         &
                  rdxsq*(f(j+1,i,k)/mapf(j+1,i)+f(j-1,i,k)/mapf(j-1,i) + &
                         f(j,i+1,k)/mapf(j,i+1)+f(j,i-1,k)/mapf(j,i-1) - &
                  d_four*f(j,i,k)/mapf(j,i))*press(j,i)
          end if
        end do
      end do
    end if
    if ( ma%hastop ) then
      i = idi2
      do k = 1 , kz
        do j = jdi1 , jdi2
          if ( ind == 0 ) then
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *         &
                  rdxsq*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k) - &
                  d_four*f(j,i,k))*press(j,i)
          else
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *         &
                  rdxsq*(f(j+1,i,k)/mapf(j+1,i)+f(j-1,i,k)/mapf(j-1,i) + &
                         f(j,i+1,k)/mapf(j,i+1)+f(j,i-1,k)/mapf(j,i-1) - &
                  d_four*f(j,i,k)/mapf(j,i))*press(j,i)
          end if
        end do
      end do
    end if
!
    call time_end(subroutine_name,idindx) 

  end subroutine diffu_d
!
  subroutine diffu_x3d(ften,f,press,xkc,kmax)
!
    implicit none
!
    integer , intent(in) :: kmax
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: xkc
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: f
    real(dp) , pointer , dimension(:,:,:) , intent(out) :: ften
    real(dp) , pointer , dimension(:,:) , intent(in) :: press
!
    integer :: i , j , k
!
    character (len=64) :: subroutine_name='diffu_x3d'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
    !
    ! fourth-order scheme for interior:
    !
    do k = 1 , kmax
      do i = icii1 , icii2
        do j = jcii1 , jcii2
          ften(j,i,k) = ften(j,i,k) - xkc(j,i,k) *    &
                      rdxsq*(f(j+2,i,k)+f(j-2,i,k) +  &
                             f(j,i+2,k)+f(j,i-2,k) -  &
                     d_four*(f(j+1,i,k)+f(j-1,i,k) +  &
                             f(j,i+1,k)+f(j,i-1,k)) + &
                    d_twelve*f(j,i,k))*press(j,i)
        end do
      end do
    end do
    !
    ! second-order scheme for east and west boundaries:
    !
    if ( ma%hasleft ) then
      j = jci1
      do k = 1 , kmax
        do i = ici1 , ici2
          ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *    &
                        rdxsq*(f(j+1,i,k)+f(j-1,i,k) + &
                               f(j,i+1,k)+f(j,i-1,k) - &
                        d_four*f(j,i,k))*press(j,i)
        end do
      end do
    end if
    if ( ma%hasright ) then
      j = jci2
      do k = 1 , kmax
        do i = ici1 , ici2
          ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *  &
                        rdxsq*(f(j+1,i,k)+f(j-1,i,k) + &
                               f(j,i+1,k)+f(j,i-1,k) - &
                        d_four*f(j,i,k))*press(j,i)
        end do
      end do
    end if
    !
    ! second-order scheme for north and south boundaries:
    !
    if ( ma%hasbottom ) then
      i = ici1
      do k = 1 , kmax
        do j = jci1 , jci2
          ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *  &
                      rdxsq*(f(j+1,i,k)+f(j-1,i,k) + &
                             f(j,i+1,k)+f(j,i-1,k) - &
                      d_four*f(j,i,k))*press(j,i)
        end do
      end do
    end if
    if ( ma%hastop ) then
      i = ici2
      do k = 1 , kmax
        do j = jci1 , jci2
          ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) *  &
                      rdxsq*(f(j+1,i,k)+f(j-1,i,k) + &
                             f(j,i+1,k)+f(j,i-1,k) - &
                      d_four*f(j,i,k))*press(j,i)
        end do
      end do
    end if

    call time_end(subroutine_name,idindx)

  end subroutine diffu_x3d
!
  subroutine diffu_x4d(ften,f,press,xkc,kmax)
!
    implicit none
!
    integer , intent(in) :: kmax
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: xkc
    real(dp) , pointer , dimension(:,:,:,:) , intent(in) :: f
    real(dp) , pointer , dimension(:,:,:,:) , intent(out) :: ften
    real(dp) , pointer , dimension(:,:) , intent(in) :: press
!
    integer :: i , j , k , n
!
    character (len=64) :: subroutine_name='diffu_x4d'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
    !
    ! fourth-order scheme for interior:
    !
    do n = 1 , ntr
      do k = 1 , kmax
        do i = icii1 , icii2
          do j = jcii1 , jcii2
            ften(j,i,k,n) = ften(j,i,k,n) - xkc(j,i,k) *    &
                        rdxsq*(f(j+2,i,k,n)+f(j-2,i,k,n) +  &
                               f(j,i+2,k,n)+f(j,i-2,k,n) -  &
                       d_four*(f(j+1,i,k,n)+f(j-1,i,k,n) +  &
                               f(j,i+1,k,n)+f(j,i-1,k,n)) + &
                      d_twelve*f(j,i,k,n))*press(j,i)
          end do
        end do
      end do
    end do
    !
    ! second-order scheme for east and west boundaries:
    !
    if ( ma%hasleft ) then
      j = jci1
      do n = 1 , ntr
        do k = 1 , kmax
          do i = ici1 , ici2
            ften(j,i,k,n) = ften(j,i,k,n) + xkc(j,i,k) *     &
                          rdxsq*(f(j+1,i,k,n)+f(j-1,i,k,n) + &
                                 f(j,i+1,k,n)+f(j,i-1,k,n) - &
                          d_four*f(j,i,k,n))*press(j,i)
          end do
        end do
      end do
    end if
    if ( ma%hasright ) then
      j = jci2
      do n = 1 , ntr
        do k = 1 , kmax
          do i = ici1 , ici2
            ften(j,i,k,n) = ften(j,i,k,n) + xkc(j,i,k) *     &
                          rdxsq*(f(j+1,i,k,n)+f(j-1,i,k,n) + &
                                 f(j,i+1,k,n)+f(j,i-1,k,n) - &
                          d_four*f(j,i,k,n))*press(j,i)
          end do
        end do
      end do
    end if
    !
    ! second-order scheme for north and south boundaries:
    !
    if ( ma%hasbottom ) then
      i = ici1
      do n = 1 , ntr
        do k = 1 , kmax
          do j = jci1 , jci2
            ften(j,i,k,n) = ften(j,i,k,n) + xkc(j,i,k) *   &
                        rdxsq*(f(j+1,i,k,n)+f(j-1,i,k,n) + &
                               f(j,i+1,k,n)+f(j,i-1,k,n) - &
                        d_four*f(j,i,k,n))*press(j,i)
          end do
        end do
      end do
    end if
    if ( ma%hastop ) then
      i = ici2
      do n = 1 , ntr
        do k = 1 , kmax
          do j = jci1 , jci2
            ften(j,i,k,n) = ften(j,i,k,n) + xkc(j,i,k) *   &
                        rdxsq*(f(j+1,i,k,n)+f(j-1,i,k,n) + &
                               f(j,i+1,k,n)+f(j,i-1,k,n) - &
                        d_four*f(j,i,k,n))*press(j,i)
          end do
        end do
      end do
    end if
    call time_end(subroutine_name,idindx)
  end subroutine diffu_x4d

end module mod_diffusion
