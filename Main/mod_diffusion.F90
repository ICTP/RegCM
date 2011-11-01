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
  use mod_runparams
  use mod_service 
  private

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
!     j       : j'th slice of variable bd3d on dot points             c
!                                                                     c
!     ind = 1 : var is already multiplied by map scale factor         c
!         = 0 : var is "not"   multiplied by map scale factor         c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine diffu_d(jstart,jend,istart,iend,ften,bd3d,press,mapf,xkc,ind)
!
    implicit none
!
    integer , intent(in) :: ind , jstart , jend , istart , iend
    real(8) , pointer , dimension(:,:,:) , intent(in) :: bd3d
    real(8) , pointer , dimension(:,:,:) , intent(in) :: xkc
    real(8) , pointer , dimension(:,:) , intent(in) :: press
    real(8) , pointer , dimension(:,:) , intent(in) :: mapf
    real(8) , pointer , dimension(:,:,:) , intent(out) :: ften
!
    integer :: i , j , k , jm1 , jm2 , jp1, jp2
!
    character (len=64) :: subroutine_name='diffu_d'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)

    do j = jstart , jend
      jm1 = j - 1
      jm2 = j - 2
      jp1 = j + 1
      jp2 = j + 2

#ifdef BAND
!
!---------------------------------------------------------------------
!
!   fourth-order scheme for interior:
!
      do k = 1 , kz
        do i = istart+1 , iend-1
          if ( ind == 0 ) then
            ften(i,k,j) = ften(i,k,j) - xkc(i,k,j) *          &
                        rdxsq*(bd3d(i,k,jp2)+bd3d(i,k,jm2) +  &
                        bd3d(i+2,k,j)+bd3d(i-2,k,j) -         &
                        d_four*(bd3d(i,k,jp1)+bd3d(i,k,jm1) + &
                        bd3d(i+1,k,j)+bd3d(i-1,k,j)) +        &
                        d_twelve*bd3d(i,k,j))*press(i,j)
          else
            ften(i,k,j) = ften(i,k,j) - xkc(i,k,j) *        &
                        rdxsq*(bd3d(i,k,jp2)/mapf(i,jp2) +  &
                        bd3d(i,k,jm2)/mapf(i,jm2) +         &
                        bd3d(i+2,k,j)/mapf(i+2,j) +         &
                        bd3d(i-2,k,j)/mapf(i-2,j) -         &
                        d_four*(bd3d(i,k,jp1)/mapf(i,jp1) + &
                        bd3d(i,k,jm1)/mapf(i,jm1) +         &
                        bd3d(i+1,k,j)/mapf(i+1,j) +         &
                        bd3d(i-1,k,j)/mapf(i-1,j)) +        &
                        d_twelve*bd3d(i,k,j)/mapf(i,j))*press(i,j)
          end if
        end do
      end do
!
!   second-order scheme for north and south boundaries:
!
      do i = istart , iend , iend-2
        do k = 1 , kz
          if ( ind == 0 ) then
            ften(i,k,j) = ften(i,k,j) + xkc(i,k,j) *         &
                        rdxsq*(bd3d(i,k,jp1)+bd3d(i,k,jm1) + &
                        bd3d(i+1,k,j)+bd3d(i-1,k,j) -        &
                        d_four*bd3d(i,k,j))*press(i,j)
          else
            ften(i,k,j) = ften(i,k,j) + xkc(i,k,j) *         &
                        rdxsq*(bd3d(i,k,jp1)/mapf(i,jp1) +   &
                        bd3d(i,k,jm1)/mapf(i,jm1) +          &
                        bd3d(i+1,k,j)/mapf(i+1,j) +          &
                        bd3d(i-1,k,j)/mapf(i-1,j) -          &
                        d_four*bd3d(i,k,j)/mapf(i,j))*press(i,j)
          end if
        end do
      end do

#else
!
!---------------------------------------------------------------------
!
      if ( (myid == 0 .and. j == 2) .or. &
           (myid == nproc-1 .and. j == jendx) ) then
!
!     second-order scheme for east or west boundary:
!
        do k = 1 , kz
          do i = istart , iend
            if ( ind == 0 ) then
              ften(i,k,j) = ften(i,k,j) + xkc(i,k,j) *         &
                          rdxsq*(bd3d(i,k,jp1)+bd3d(i,k,jm1) + &
                          bd3d(i+1,k,j)+bd3d(i-1,k,j) -        &
                          d_four*bd3d(i,k,j))*press(i,j)
            else
              ften(i,k,j) = ften(i,k,j) + xkc(i,k,j) *       &
                          rdxsq*(bd3d(i,k,jp1)/mapf(i,jp1) + &
                          bd3d(i,k,jm1)/mapf(i,jm1) +        &
                          bd3d(i+1,k,j)/mapf(i+1,j) +        &
                          bd3d(i-1,k,j)/mapf(i-1,j) -        &
                          d_four*bd3d(i,k,j)/mapf(i,j))*press(i,j)
            end if
          end do
        end do
!
      else
!
!     fourth-order scheme for interior:
!
        do k = 1 , kz
          do i = istart+1 , iend-1
            if ( ind == 0 ) then
              ften(i,k,j) = ften(i,k,j) - xkc(i,k,j) *          &
                          rdxsq*(bd3d(i,k,jp2)+bd3d(i,k,jm2) +  &
                          bd3d(i+2,k,j)+bd3d(i-2,k,j) -         &
                          d_four*(bd3d(i,k,jp1)+bd3d(i,k,jm1) + &
                          bd3d(i+1,k,j)+bd3d(i-1,k,j)) +        &
                          d_twelve*bd3d(i,k,j))*press(i,j)
            else
              ften(i,k,j) = ften(i,k,j) - xkc(i,k,j) *        &
                          rdxsq*(bd3d(i,k,jp2)/mapf(i,jp2) +  &
                          bd3d(i,k,jm2)/mapf(i,jm2) +         &
                          bd3d(i+2,k,j)/mapf(i+2,j) +         &
                          bd3d(i-2,k,j)/mapf(i-2,j) -         &
                          d_four*(bd3d(i,k,jp1)/mapf(i,jp1) + &
                          bd3d(i,k,jm1)/mapf(i,jm1) +         &
                          bd3d(i+1,k,j)/mapf(i+1,j) +         &
                          bd3d(i-1,k,j)/mapf(i-1,j)) +        &
                          d_twelve*bd3d(i,k,j)/mapf(i,j))*press(i,j)
            end if
          end do
        end do
!
!     second-order scheme for north and south boundaries:
!
        do i = istart , iend , iend-2
          do k = 1 , kz
            if ( ind == 0 ) then
              ften(i,k,j) = ften(i,k,j) + xkc(i,k,j) *         &
                          rdxsq*(bd3d(i,k,jp1)+bd3d(i,k,jm1) + &
                          bd3d(i+1,k,j)+bd3d(i-1,k,j) -        &
                          d_four*bd3d(i,k,j))*press(i,j)
            else
              ften(i,k,j) = ften(i,k,j) + xkc(i,k,j) *       &
                          rdxsq*(bd3d(i,k,jp1)/mapf(i,jp1) + &
                          bd3d(i,k,jm1)/mapf(i,jm1) +        &
                          bd3d(i+1,k,j)/mapf(i+1,j) +        &
                          bd3d(i-1,k,j)/mapf(i-1,j) -        &
                          d_four*bd3d(i,k,j)/mapf(i,j))*press(i,j)
            end if
          end do
        end do
!
      end if
#endif
!
    end do
    call time_end(subroutine_name,idindx) 

  end subroutine diffu_d
!
  subroutine diffu_x(jstart,jend,istart,iend,ften,bc3d,press,xkc,kmax)
!
    implicit none
!
    integer , intent(in) :: jstart ,jend , istart , iend , kmax
    real(8) , pointer , dimension(:,:,:) , intent(in) :: xkc
    real(8) , pointer , dimension(:,:,:) , intent(in) :: bc3d
    real(8) , pointer , dimension(:,:,:) , intent(out) :: ften
    real(8) , pointer , dimension(:,:) , intent(in) :: press
!
    integer :: i , j , k , jm1 , jm2 , jp1, jp2
!
    character (len=64) :: subroutine_name='diffu_x'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!   compute the diffusion term:
!
    do j = jstart , jend
      jm1 = j - 1
      jm2 = j - 2
      jp1 = j + 1
      jp2 = j + 2

#ifdef BAND
!
!---------------------------------------------------------------------
!
!   fourth-order scheme for interior:
!
      do k = 1 , kmax
        do i = istart+1 , iend-1
          ften(i,k,j) = ften(i,k,j) - xkc(i,k,j) *           &
                      rdxsq*(bc3d(i,k,jp2)+bc3d(i,k,jm2) +   &
                             bc3d(i+2,k,j)+bc3d(i-2,k,j) -   &
                      d_four*(bc3d(i,k,jp1)+bc3d(i,k,jm1) +  &
                              bc3d(i+1,k,j)+bc3d(i-1,k,j)) + &
                      d_twelve*bc3d(i,k,j))*press(i,j)
        end do
      end do
!
!   second-order scheme for north and south boundaries:
!
      i = istart
      do k = 1 , kmax
        ften(i,k,j) = ften(i,k,j) + xkc(i,k,j) *         &
                    rdxsq*(bc3d(i,k,jp1)+bc3d(i,k,jm1) + &
                           bc3d(i+1,k,j)+bc3d(i-1,k,j) - &
                    d_four*bc3d(i,k,j))*press(i,j)
      end do
      i = iend
      do k = 1 , kmax
        ften(i,k,j) = ften(i,k,j) + xkc(i,k,j) *         &
                    rdxsq*(bc3d(i,k,jp1)+bc3d(i,k,jm1) + &
                           bc3d(i+1,k,j)+bc3d(i-1,k,j) - &
                    d_four*bc3d(i,k,j))*press(i,j)
      end do

#else

!
!----------------------------------------------------------------------
!
      if ( (myid == 0 .and. j == 2) .or. &
           (myid == nproc-1 .and. j == jendm) ) then
!
!     second-order scheme for east or west boundary:
!
        do k = 1 , kmax
          do i = istart , iend
            ften(i,k,j) = ften(i,k,j) + xkc(i,k,j) *         &
                        rdxsq*(bc3d(i,k,jp1)+bc3d(i,k,jm1) + &
                               bc3d(i+1,k,j)+bc3d(i-1,k,j) - &
                        d_four*bc3d(i,k,j))*press(i,j)
          end do
        end do
!
      else
!
!     fourth-order scheme for interior:
!
        do k = 1 , kmax
          do i = istart+1 , iend-1
            ften(i,k,j) = ften(i,k,j) - xkc(i,k,j) *           &
                        rdxsq*(bc3d(i,k,jp2)+bc3d(i,k,jm2) +   &
                               bc3d(i+2,k,j)+bc3d(i-2,k,j) -   &
                        d_four*(bc3d(i,k,jp1)+bc3d(i,k,jm1) +  &
                                bc3d(i+1,k,j)+bc3d(i-1,k,j)) + &
                        d_twelve*bc3d(i,k,j))*press(i,j)
          end do
        end do
!
!     second-order scheme for north and south boundaries:
!
        i = istart
        do k = 1 , kmax
          ften(i,k,j) = ften(i,k,j) + xkc(i,k,j) *         &
                      rdxsq*(bc3d(i,k,jp1)+bc3d(i,k,jm1) + &
                             bc3d(i+1,k,j)+bc3d(i-1,k,j) - &
                      d_four*bc3d(i,k,j))*press(i,j)
        end do
        i = iend
        do k = 1 , kmax
          ften(i,k,j) = ften(i,k,j) + xkc(i,k,j) *         &
                      rdxsq*(bc3d(i,k,jp1)+bc3d(i,k,jm1) + &
                             bc3d(i+1,k,j)+bc3d(i-1,k,j) - &
                      d_four*bc3d(i,k,j))*press(i,j)
        end do
!
      end if
!
#endif
    end do

    call time_end(subroutine_name,idindx)

  end subroutine diffu_x
!
end module mod_diffusion
