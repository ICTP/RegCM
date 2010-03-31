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

      subroutine smther(slab,is1,is2,npass,point,iflg)
      implicit none
!
! Dummy arguments
!
      integer :: iflg , is1 , is2 , npass
      character(5) :: point
      real(4) , dimension(is1,is2) :: slab
      intent (in) iflg , is1 , is2 , npass , point
      intent (inout) slab
!
! Local variables
!
      real(4) :: aplus , asv , cell
      integer :: i , icross , ie , iem , j , je , jem , k , kp
      real(4) , dimension(2) :: xnu
!
!     purpose: spatially smooth data in slab to dampen short
!     wavelength components
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
!         first smooth in the is1 direction
          do i = 1 , ie
            asv = slab(i,1)
            do j = 2 , je - 1
              aplus = slab(i,j+1)
              cell = slab(i,j)
              slab(i,j) = slab(i,j) + xnu(kp)                           &
                        & *((asv+aplus)/2.0-slab(i,j))
              if ( iflg==0 ) then
                if ( (i>6) .and. (i<iem) .and. (j>6) .and. (j<jem) )    &
                   & slab(i,j) = cell
              else if ( iflg==1 ) then
                if ( i>20 ) slab(i,j) = cell
              else
              end if
              asv = cell
            end do
          end do
!         smooth in the is2 direction
          do j = 1 , je
            asv = slab(1,j)
            do i = 2 , ie - 1
              aplus = slab(i+1,j)
              cell = slab(i,j)
              slab(i,j) = slab(i,j) + xnu(kp)                           &
                        & *((asv+aplus)/2.0-slab(i,j))
              if ( iflg==0 ) then
                if ( (i>6) .and. (i<iem) .and. (j>6) .and. (j<jem) )    &
                   & slab(i,j) = cell
              else if ( iflg==1 ) then
                if ( i>20 ) slab(i,j) = cell
              else
              end if
              asv = cell
            end do
          end do
!         40 continue
        end do
      end do
      end subroutine smther
