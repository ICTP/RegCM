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

module mod_smooth

  use mod_intkinds
  use mod_realkinds
  use mod_constants

  contains

  subroutine smth121(htgrid,jx,iy)
    implicit none
    integer(ik4) , intent(in) :: jx , iy
    real(rk8) , intent(inout) , dimension(jx,iy) :: htgrid
!
    integer(ik4) :: i , j
    real(rk8) , dimension(jx,iy) :: hscr1
    !
    ! PURPOSE :  PERFORMS THE 1-2-1 SMOOTHING TO REMOVE PRIMARILY THE
    ! 2DX WAVES FROM THE FIELDS htgrid
    !
    hscr1(:,:) = htgrid(:,:)
    do i = 1 , iy
      do j = 2 , jx - 1
        if ( (htgrid(j,i) <= -d_one) .or. (htgrid(j,i) > d_zero) ) then
          hscr1(j,i) = d_rfour*(d_two*htgrid(j,i)+htgrid(j+1,i)+htgrid(j-1,i))
        end if
      end do
    end do
    do i = 2 , iy - 1
      do j = 1 , jx
        if ( (hscr1(j,i) <= -d_one) .or. (hscr1(j,i) > d_zero) ) then
          htgrid(j,i) = d_rfour*(d_two*hscr1(j,i)+hscr1(j,i+1)+hscr1(j,i-1))
        end if
      end do
    end do
  end subroutine smth121

  subroutine smther(slab,nj,ni,npass,iflg)
    implicit none
    integer(ik4) , intent(in) :: iflg , ni , nj , npass
    real(rk8) , intent(inout) , dimension(nj,ni) :: slab
!
    real(rk8) :: aplus , asv , cell
    integer(ik4) :: i , ie , iem , j , je , jem , k , kp
    real(rk8) , dimension(2) :: xnu
    !
    ! purpose: spatially smooth data in slab to dampen short
    ! wavelength components
    !
    ie = ni
    je = nj
    iem = ie - 5
    jem = je - 5
    xnu(1) =  d_half
    xnu(2) = -d_half
    do k = 1 , npass
      do kp = 1 , 2
        ! first smooth in the ni direction
        do i = 1 , ie - 1
          asv = slab(1,i)
          do j = 2 , je - 1
            aplus = slab(j,i+1)
            cell = slab(j,i)
            slab(j,i) = cell + xnu(kp)*( (asv+aplus)/d_two - cell)
            if ( iflg == 0 ) then
              if ( (i > 6) .and. (i < iem) .and. &
                   (j > 6) .and. (j < jem) ) then
                slab(j,i) = cell
              end if
            else if ( iflg == 1 ) then
              if ( i > 20 ) slab(j,i) = cell
            end if
            asv = cell
          end do
        end do
        ! smooth in the nj direction
        do j = 1 , je - 1
          asv = slab(j,1)
          do i = 2 , ie - 1
            aplus = slab(j,i+1)
            cell = slab(j,i)
            slab(j,i) = cell + xnu(kp)*((asv+aplus)/d_two - cell)
            if ( iflg == 0 ) then
              if ( (i > 6) .and. (i < iem) .and. &
                   (j > 6) .and. (j < jem) ) then
                slab(j,i) = cell
              end if
            else if ( iflg == 1 ) then
              if ( i > 20 ) slab(j,i) = cell
            end if
            asv = cell
          end do
        end do
      end do
    end do
  end subroutine smther

  subroutine smthtr(slab1,nj,ni)
    implicit none
    integer(ik4) , intent(in) :: nj , ni
    real(rk8) , intent(inout) , dimension(nj,ni) :: slab1
!
    integer(ik4) :: i , iflg , j , k , n , n1 , npass
    integer(ik4) , allocatable , dimension(:) :: ii , jj
    !
    ! smooth terrain arrays
    !
    allocate(ii(nj*ni))
    allocate(jj(nj*ni))
    n = 1
    do i = 1 , ni
      do j = 1 , nj
        if ( slab1(j,i) < d_zero ) then
          ii(n) = i
          jj(n) = j
          slab1(j,i) = d_zero
          n = n + 1
        end if
      end do
    end do
    n1 = n - 1
    npass = 10
    iflg = 0          ! 0 = smoothing only at boundary
    call smther(slab1,nj,ni,npass,iflg)
!   npass = 2000
!   iflg = 0          ! 1 = extensive smoothing over south boundary
!   call smther(slab1,nj,ni,npass,iflg)
    do i = 1 , ni
      do j = 1 , nj
        slab1(j,i) = slab1(j,i)
        if ( slab1(j,i) < d_zero ) slab1(j,i) = d_zero
      end do
    end do
    do k = 1 , n - 1
      i = ii(k)
      j = jj(k)
      slab1(j,i) = -0.001D0
    end do
    deallocate(ii)
    deallocate(jj)
  end subroutine smthtr

end module mod_smooth
