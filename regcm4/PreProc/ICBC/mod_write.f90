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

      module mod_write

      implicit none

      integer :: noutrec
      integer :: iny , jnx , knz

      real(4) , allocatable , dimension(:,:) :: ps4 , ts4
      real(4) , allocatable , dimension(:,:,:) :: c4 , h4 , q4
      real(4) , allocatable , dimension(:,:,:) :: t4 , u4 , v4

      data noutrec /0/

      contains

      subroutine init_output(jx,iy,kz)
      implicit none
      integer , intent(in) :: jx , iy , kz
      iny = iy
      jnx = jx
      knz = kz
      allocate(ps4(jnx,iny))
      allocate(ts4(jnx,iny))
      allocate(c4(jnx,iny,knz))
      allocate(h4(jnx,iny,knz))
      allocate(q4(jnx,iny,knz))
      allocate(t4(jnx,iny,knz))
      allocate(u4(jnx,iny,knz))
      allocate(v4(jnx,iny,knz))
      end subroutine init_output

      subroutine free_output
      implicit none
      deallocate(ps4)
      deallocate(ts4)
      deallocate(c4)
      deallocate(h4)
      deallocate(q4)
      deallocate(t4)
      deallocate(u4)
      deallocate(v4)
      end subroutine free_output

      subroutine writef(ptop,idate)
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
      write (64,rec=noutrec) idate , jnx , iny , knz
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((u4(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((v4(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((t4(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((q4(i,j,k),i=1,jnx),j=1,iny)
      end do
      noutrec = noutrec + 1
      write (64,rec=noutrec) ((ps4(i,j)+ptop,i=1,jnx),j=1,iny)
      noutrec = noutrec + 1
      write (64,rec=noutrec) ts4
!
      end subroutine writef

      subroutine writef2(ptop,idate)
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
      write (64,rec=noutrec) idate , jnx , iny , knz
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((u4(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((v4(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((t4(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((q4(i,j,k),i=1,jnx),j=1,iny)
      end do
      noutrec = noutrec + 1
      write (64,rec=noutrec) ((ps4(i,j)+ptop,i=1,jnx),j=1,iny)
      noutrec = noutrec + 1
      write (64,rec=noutrec) ts4
!
      end subroutine writef2

      subroutine writefs(qs3,ti3,ts3,snow,ptop,idate,lsmtyp)
      implicit none
!
! Dummy arguments
!
      integer :: idate
      character(4) :: lsmtyp
      real(4) :: ptop
      real(4) , dimension(jnx,iny,4) :: qs3 , ti3 , ts3
      real(4) , dimension(jnx,iny) :: snow
      intent (in) idate , lsmtyp , ptop , snow , ti3 , ts3 , qs3
!
! Local variables
!
      integer :: i , j , k
!
!     THIS ROUTINE WRITES OUT AN INPUT FILE FOR THE RCM
!
      noutrec = noutrec + 1
      write (64,rec=noutrec) idate , jnx , iny , knz
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((u4(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((v4(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((t4(i,j,k),i=1,jnx),j=1,iny)
      end do
      do k = knz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((q4(i,j,k),i=1,jnx),j=1,iny)
      end do
      noutrec = noutrec + 1
      write (64,rec=noutrec) ((ps4(i,j)+ptop,i=1,jnx),j=1,iny)
      noutrec = noutrec + 1
      write (64,rec=noutrec) ts4
 
      if ( lsmtyp=='USGS' ) then
        do k = 1 , 4
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((qs3(i,j,k),i=1,jnx),j=1,iny)
        end do
        do k = 1 , 4
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((ti3(i,j,k),i=1,jnx),j=1,iny)
        end do
        do k = 1 , 4
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((ts3(i,j,k),i=1,jnx),j=1,iny)
        end do
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((snow(i,j),i=1,jnx),j=1,iny)
      end if
!
      end subroutine writefs

      end module mod_write
