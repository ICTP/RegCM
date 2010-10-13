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

      module mod_wrtoxd

      implicit none
      integer :: noutrec
      integer :: iny , jnx , knz
      real(4) , allocatable , dimension(:,:,:) :: oh4,ho24,o34,no34,h2o24

      data noutrec /0/

      contains

      subroutine init_outoxd(jx,iy,kz)
      implicit none
      integer , intent(in) :: jx , iy , kz
      iny = iy
      jnx = jx
      knz = kz
      allocate(oh4(jnx,iny,knz))
      allocate(ho24(jnx,iny,knz))
      allocate(o34(jnx,iny,knz))
      allocate(no34(jnx,iny,knz))
      allocate(h2o24(jnx,iny,knz))
      end subroutine init_outoxd

      subroutine free_outoxd
      implicit none
      deallocate(oh4)
      deallocate(ho24)
      deallocate(o34)
      deallocate(no34)
      deallocate(h2o24)
      end subroutine free_outoxd

      subroutine writeox(ptop,idate)
      implicit none
!
! Dummy arguments
!
      integer :: idate
      real(4) :: ptop
      intent (in) idate , ptop
!
! Local variables
!
      integer :: i , j , k
!
!     THIS ROUTINE WRITES OUT AN INPUT FILE FOR THE RCM
!
!     PRINT *,'WRITING OUTPUT:  IDATE= ',IDATE

      noutrec = noutrec + 1
      write(*,*)'INSID',noutrec,idate
      write (64,rec=noutrec) idate , jnx , iny , knz
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((oh4(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((ho24(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((o34(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((no34(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((h2o24(i,j,k),i=1,jnx),j=1,iny)
      end do
!
      end subroutine writeox

      end module mod_wrtoxd
