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

  public :: smth121 , smthtr , cleanpeaks , h2oelev , remove_high_gradients

  contains

  subroutine smth121(htgrid,jx,iy)
    implicit none
    integer(ik4) , intent(in) :: jx , iy
    real(rkx) , intent(inout) , dimension(jx,iy) :: htgrid
    integer(ik4) :: i , j
    real(rkx) , dimension(jx,iy) :: hscr1 , hscr2
    !
    ! PURPOSE :  PERFORMS THE 1-2-1 SMOOTHING TO REMOVE PRIMARILY THE
    ! 2DX WAVES FROM THE FIELDS htgrid.
    !
    hscr1(:,:) = htgrid(:,:)
    hscr2(:,:) = htgrid(:,:)
    do i = 2 , iy - 1
      do j = 2 , jx - 1
        hscr2(j,i) = d_rfour*(d_two*hscr1(j,i)+hscr1(j+1,i)+hscr1(j-1,i))
      end do
    end do
    do i = 2 , iy - 1
      do j = 2 , jx - 1
        hscr1(j,i) = d_rfour*(d_two*hscr2(j,i)+hscr2(j,i+1)+hscr2(j,i-1))
      end do
    end do
    htgrid(:,:) = hscr1(:,:)
  end subroutine smth121

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
    ie = ni
    je = nj
    iem = ie - nsp + 1
    jem = je - nsp + 1
    xnu(1) =  d_half
    xnu(2) = -d_half
    do k = 1 , npass
      do kp = 1 , 2
        ! first smooth in the ni direction
        do i = 1 , ie - 1
          asv = slab(1,i)
          do j = 2 , je - 1
            if ( (i <= nsp) .and. (i >= iem) .and. &
                 (j <= nsp) .and. (j >= jem) ) then
              cell = slab(j,i)
              aplus = slab(j,i+1)
              slab(j,i) = cell + xnu(kp)*( (asv+aplus)/d_two - cell)
            end if
            asv = cell
          end do
        end do
        ! smooth in the nj direction
        do j = 1 , je - 1
          asv = slab(j,1)
          do i = 2 , ie - 1
            if ( (i <= nsp) .and. (i >= iem) .and. &
                 (j <= nsp) .and. (j >= jem) ) then
              cell = slab(j,i)
              aplus = slab(j,i+1)
              slab(j,i) = cell + xnu(kp)*((asv+aplus)/d_two - cell)
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
    npass = 10
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

  subroutine h2oelev(nj,ni,xt,mask)
    implicit none
    integer(ik4) , intent(in) :: nj , ni
    real(rkx) , intent(in) , dimension(nj,ni) :: mask
    real(rkx) , intent(inout) , dimension(nj,ni) :: xt
    real(rkx) , dimension(nj,ni) :: xtn
    real(rkx) :: np , sump
    integer(ik4) :: i , j , ii , jj , iii , jjj , ipass
    integer , parameter :: npass = 1

    do ipass = 1 , npass
      do i = 1 , ni
        do j = 1 , nj
          if ( mask(j,i) < 1.0_rkx ) then
            np = 4.0_rkx
            sump = 4.0_rkx * xt(j,i)
            do ii = i-1 , i+1
              do jj = j-1 , j+1
                iii = max(1,min(ii,ni))
                jjj = max(1,min(jj,nj))
                sump = sump + xt(jjj,iii)
                np = np + 1
              end do
            end do
            xtn(j,i) = sump/np
          else
            xtn(j,i) = xt(j,i)
          end if
        end do
      end do
      xt(:,:) = xtn(:,:)
    end do
  end subroutine h2oelev

  subroutine remove_high_gradients(nj,ni,xt)
    implicit none
    integer(ik4) , intent(in) :: nj , ni
    real(rkx) , intent(inout) , dimension(nj,ni) :: xt
    real(rkx) , dimension(nj,ni) :: xtn
    real(rkx) , dimension(nj,ni) :: dhdx1
    real(rkx) , dimension(nj,ni) :: dhdx2
    real(rkx) , dimension(nj,ni) :: dhdy1
    real(rkx) , dimension(nj,ni) :: dhdy2
    real(rkx) :: np , sump
    integer(ik4) :: i , j , ii , jj , iii , jjj , ipass
    integer , parameter :: npass = 2
    real(rkx) , parameter :: maxgra = 750.0_rkx
    dhdx1(:,:) = d_zero
    dhdy1(:,:) = d_zero
    dhdx2(:,:) = d_zero
    dhdy2(:,:) = d_zero
    do i = 2 , ni-1
      do j = 2 , nj-1
        dhdx1(j,i) = abs(xt(j,i) - xt(j-1,i))
        dhdx2(j,i) = abs(xt(j+1,i) - xt(j,i))
        dhdy1(j,i) = abs(xt(j,i) - xt(j,i-1))
        dhdy2(j,i) = abs(xt(j,i+1) - xt(j,i))
      end do
    end do
    do ipass = 1 , npass
      do i = 1 , ni
        do j = 1 , nj
          if ( dhdx1(j,i) > maxgra .or. &
               dhdx2(j,i) > maxgra .or. &
               dhdy1(j,i) > maxgra .or. &
               dhdy2(j,i) > maxgra ) then
            np = 0.0_rkx
            sump = 0.0_rkx
            do ii = i-1 , i+1
              do jj = j-1 , j+1
                iii = max(1,min(ii,ni))
                jjj = max(1,min(jj,nj))
                np = np + 1
                sump = sump + xt(jjj,iii)
              end do
            end do
            if ( np > 0.0_rkx ) then
              xtn(j,i) = xt(j,i) + (sump - xt(j,i) * np)/np
            else
              xtn(j,i) = xt(j,i)
            end if
          else
            xtn(j,i) = xt(j,i)
          end if
        end do
      end do
      xt(:,:) = xtn(:,:)
    end do
  end subroutine remove_high_gradients

  subroutine cleanpeaks(nj,ni,xt)
    implicit none
    integer(ik4) , intent(in) :: nj , ni
    real(rkx) , intent(inout) , dimension(nj,ni) :: xt
    real(rkx) , dimension(nj,ni) :: xtn
    real(rkx) :: np , sump
    integer(ik4) :: i , j , ii , jj , iii , jjj , ipass
    integer , parameter :: npass = 1
    real(rkx) , parameter :: mhgt = 5.0_rkx

    do ipass = 1 , npass
      do i = 1 , ni
        do j = 1 , nj
          np = 0.0_rkx
          sump = 0.0_rkx
          do ii = i-1 , i+1
            do jj = j-1 , j+1
              iii = max(1,min(ii,ni))
              jjj = max(1,min(jj,nj))
              if ( xt(jjj,iii) > mhgt ) then
                np = np + 1
                sump = sump + xt(jjj,iii)
              end if
            end do
          end do
          if ( np > 0.0_rkx ) then
            xtn(j,i) = xt(j,i) + (sump - xt(j,i) * np)/np
          else
            xtn(j,i) = xt(j,i)
          end if
        end do
      end do
      xt(:,:) = xtn(:,:)
    end do
  end subroutine cleanpeaks

end module mod_smooth
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
