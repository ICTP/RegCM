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
  use mod_stdio

  contains

  subroutine smth121(htgrid,iy,jx)
    implicit none
!
    integer(ik4) , intent(in) :: iy , jx
    real(rk4) , intent(inout) , dimension(iy,jx) :: htgrid
!
    integer(ik4) :: i , j
    real(rk4) , dimension(iy,jx) :: hscr1
    !
    ! PURPOSE :  PERFORMS THE 1-2-1 SMOOTHING TO REMOVE PRIMARILY THE
    ! 2DX WAVES FROM THE FIELDS htgrid
    !
    do j = 1 , jx
      do i = 1 , iy
        hscr1(i,j) = htgrid(i,j)
      end do
    end do
    do i = 1 , iy
      do j = 2 , jx - 1
        if ( (htgrid(i,j)<=-.1) .or. (htgrid(i,j)>0.) ) hscr1(i,j)    &
             = .25*(2.*htgrid(i,j)+htgrid(i,j+1)+htgrid(i,j-1))
      end do
    end do
    do j = 1 , jx
      do i = 2 , iy - 1
        if ( (hscr1(i,j)<=-.1) .or. (hscr1(i,j)>0.) ) htgrid(i,j)     &
             = .25*(2.*hscr1(i,j)+hscr1(i+1,j)+hscr1(i-1,j))
      end do
    end do
  end subroutine smth121

  subroutine smther(slab,is1,is2,npass,point,iflg)
    implicit none
!
    integer(ik4) , intent(in) :: iflg , is1 , is2 , npass
    character(5) , intent(in) :: point
    real(rk4) , intent(inout) , dimension(is1,is2) :: slab
!
    real(rk4) :: aplus , asv , cell
    integer(ik4) :: i , icross , ie , iem , j , je , jem , k , kp
    real(rk4) , dimension(2) :: xnu
    !
    ! purpose: spatially smooth data in slab to dampen short
    ! wavelength components
    !
    icross = 0
    if ( point=='cross' ) icross = 1
    ie = is1
    je = is2
    iem = ie - 5
    jem = je - 5
    xnu(1) = 0.5
    xnu(2) = -.50
    do k = 1 , npass
      do kp = 1 , 2
        ! first smooth in the is1 direction
        do i = 1 , ie
          asv = slab(i,1)
          do j = 2 , je - 1
            aplus = slab(i,j+1)
            cell = slab(i,j)
            slab(i,j) = slab(i,j) + xnu(kp)                           &
                        *((asv+aplus)/2.0-slab(i,j))
            if ( iflg==0 ) then
              if ( (i>6) .and. (i<iem) .and. (j>6) .and. (j<jem) )    &
                   slab(i,j) = cell
            else if ( iflg==1 ) then
              if ( i>20 ) slab(i,j) = cell
            else
            end if
            asv = cell
          end do
        end do
        ! smooth in the is2 direction
        do j = 1 , je
          asv = slab(1,j)
          do i = 2 , ie - 1
            aplus = slab(i+1,j)
            cell = slab(i,j)
            slab(i,j) = slab(i,j) + xnu(kp)                           &
                        *((asv+aplus)/2.0-slab(i,j))
            if ( iflg==0 ) then
              if ( (i>6) .and. (i<iem) .and. (j>6) .and. (j<jem) )    &
                   slab(i,j) = cell
            else if ( iflg==1 ) then
              if ( i>20 ) slab(i,j) = cell
            else
            end if
            asv = cell
          end do
        end do
      end do
    end do
  end subroutine smther

  subroutine smthtr(slab1,is1,is2)
    implicit none
!
    integer(ik4) , parameter :: nocean = 20000
!
    integer(ik4) , intent(in) :: is1 , is2
    real(rk4) , intent(inout) , dimension(is1,is2) :: slab1
!
    integer(ik4) :: i , iflg , j , k , n , n1 , npass
    integer(ik4) , dimension(nocean) :: ii , jj
    character(5) :: point
    !
    ! smooth terrain arrays
    !
    n = 1
    do i = 1 , is1
      do j = 1 , is2
        if ( slab1(i,j) < 0.0 ) then
          ii(n) = i
          jj(n) = j
          slab1(i,j) = 0.0
          n = n + 1
        end if
      end do
    end do
    n1 = n - 1
    write(stdout,99001) n1
    if ( n > nocean ) write(stdout,99002)
    point = 'cross'
    npass = 10
    iflg = 0          ! 0 = smoothing only at boundary
    call smther(slab1,is1,is2,npass,point,iflg)
!   npass = 2000
!   iflg = 0          ! 1 = extensive smoothing over south boundary
!   call smther( slab1, is1, is2, npass, point, iflg )
    do i = 1 , is1
      do j = 1 , is2
        slab1(i,j) = slab1(i,j)
        if ( slab1(i,j)<0.0 ) slab1(i,j) = 0.0
      end do
    end do
    do k = 1 , n - 1
      i = ii(k)
      j = jj(k)
      slab1(i,j) = -0.001
    end do
99001 format (5x,'  there are a total of ',i5,' points of ocean')
99002 format ('0',2x,'dimension exceeded in subr smthtr')
  end subroutine smthtr

end module mod_smooth
