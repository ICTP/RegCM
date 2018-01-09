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

  private

  public :: smtdsmt , smth121 , smthtr

  contains

  subroutine smth121(htgrid,jx,iy)
    implicit none
    integer(ik4) , intent(in) :: jx , iy
    real(rkx) , intent(inout) , dimension(jx,iy) :: htgrid
    integer(ik4) :: n , i , j
    real(rkx) , dimension(jx,iy) :: hscr
    integer(ik4) , parameter :: npass = 1
    !
    ! PURPOSE :  PERFORMS THE 1-2-1 SMOOTHING TO REMOVE PRIMARILY THE
    ! 2DX WAVES FROM THE FIELDS htgrid.
    !
    hscr(:,:) = htgrid(:,:)
    do n = 1 , npass
      do i = 1 , iy
        do j = 2 , jx - 2
          hscr(j,i) = d_rfour*(d_two*htgrid(j,i)+htgrid(j+1,i)+htgrid(j-1,i))
          if ( hscr(j,i) < d_one ) hscr(j,i) = d_zero
        end do
      end do
      do i = 2 , iy - 2
        do j = 1 , jx
          htgrid(j,i) = d_rfour*(d_two*hscr(j,i)+hscr(j,i+1)+hscr(j,i-1))
          if ( htgrid(j,i) < d_one ) htgrid(j,i) = d_zero
        end do
      end do
    end do
  end subroutine smth121

  subroutine smtdsmt(slab,nj,ni)
    implicit none
    integer(ik4) , intent(in) :: ni , nj
    real(rkx) , intent(inout) , dimension(nj,ni) :: slab
    real(rkx) :: aplus , asv , cell
    integer(ik4) :: i , ie , j , je , np , kp
    real(rkx) , dimension(2) :: xnu
    integer(ik4) , parameter :: npass = 2
    !
    ! purpose: spatially smooth data in slab to dampen short
    ! wavelength components
    !
    ie = ni-1
    je = nj-1
    xnu(1) =  0.50_rkx
    xnu(2) = -0.52_rkx
    do np = 1 , npass
      do kp = 1 , 2
        ! first smooth in the ni direction
        do i = 1 , ni
          asv = slab(1,i)
          do j = 2 , je
            cell = slab(j,i)
            aplus = slab(j+1,i)
            slab(j,i) = cell + xnu(kp)*( (asv+aplus)/d_two - cell )
            if ( slab(j,i) < d_one ) slab(j,i) = d_zero
            asv = cell
          end do
        end do
        ! smooth in the nj direction
        do j = 1 , nj
          asv = slab(j,1)
          do i = 2 , ie
            cell = slab(j,i)
            aplus = slab(j,i+1)
            slab(j,i) = cell + xnu(kp)*( (asv+aplus)/d_two - cell )
            if ( slab(j,i) < d_one ) slab(j,i) = d_zero
            asv = cell
          end do
        end do
      end do
    end do
  end subroutine smtdsmt

  subroutine smther(slab,nj,ni,nsp,npass)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nsp , npass
    real(rkx) , intent(inout) , dimension(nj,ni) :: slab
    real(rkx) :: aplus , asv , cell
    integer(ik4) :: i , ie , iem , j , je , jem , k , kp
    real(rkx) , dimension(2) :: xnu
    !
    ! purpose: spatially smooth data in slab to dampen short
    ! wavelength components
    !
    ie = ni - 1
    je = nj - 1
    iem = ie - nsp + 1
    jem = je - nsp + 1
    xnu(1) =  0.50_rkx
    xnu(2) = -0.52_rkx
    do k = 1 , npass
      do kp = 1 , 2
        ! first smooth in the ni direction
        do i = 2 , ie
          asv = slab(1,i)
          do j = 2 , je
            cell = slab(j,i)
            aplus = slab(j+1,i)
            if ( (i <= nsp) .or. (i >= iem) .and. &
                 (j <= nsp) .or. (j >= jem) ) then
              slab(j,i) = max(cell + xnu(kp)*( (asv+aplus)/d_two - cell),d_zero)
            end if
            asv = cell
          end do
        end do
        ! smooth in the nj direction
        do j = 2 , je
          asv = slab(j,1)
          do i = 2 , ie
            cell = slab(j,i)
            aplus = slab(j,i+1)
            if ( (i <= nsp) .or. (i >= iem) .and. &
                 (j <= nsp) .or. (j >= jem) ) then
              slab(j,i) = max(cell + xnu(kp)*((asv+aplus)/d_two - cell),d_zero)
            end if
            asv = cell
          end do
        end do
      end do
    end do
  end subroutine smther

  subroutine smthtr(slab1,nj,ni,nsp)
    implicit none
    integer(ik4) , intent(in) :: nj , ni , nsp
    real(rkx) , intent(inout) , dimension(nj,ni) :: slab1
    integer(ik4) :: i , j , k , n , npass
    integer(ik4) , allocatable , dimension(:) :: ii , jj
    !
    ! smooth terrain arrays
    !
    ! Get points with negative elevation
    !
    allocate(ii(nj*ni))
    allocate(jj(nj*ni))
    n = 0
    do i = 1 , ni
      do j = 1 , nj
        if ( slab1(j,i) < d_zero ) then
          n = n + 1
          ii(n) = i
          jj(n) = j
          slab1(j,i) = d_zero
        end if
      end do
    end do
    !
    ! Apply filter
    !
    npass = 2
    call smther(slab1,nj,ni,nsp,npass)
    !
    ! Remove negative elevation points
    !
    do i = 1 , ni
      do j = 1 , nj
        if ( slab1(j,i) < d_zero ) slab1(j,i) = d_zero
      end do
    end do
    !
    ! Put back points with negative elevation
    !
    do k = 1 , n
      i = ii(k)
      j = jj(k)
      slab1(j,i) = -0.001_rkx
    end do
    deallocate(ii)
    deallocate(jj)
  end subroutine smthtr

end module mod_smooth
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
