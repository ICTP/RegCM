!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_write

      use mod_regcm_param , only : jx , iy , kz

      implicit none

      integer :: noutrec

      real , dimension(jx,iy) :: ps4 , ts4
      real , dimension(jx,iy,kz) :: c4 , h4 , q4 , t4 , u4 , v4

      data noutrec /0/

      contains

      subroutine writef(ptop,idate)
      implicit none
!
! Dummy arguments
!
      integer :: idate
      real :: ptop
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
      write (64,rec=noutrec) idate , jx , iy , kz
      do k = kz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((u4(i,j,k),i=1,jx),j=1,iy)
      end do
      do k = kz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((v4(i,j,k),i=1,jx),j=1,iy)
      end do
      do k = kz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((t4(i,j,k),i=1,jx),j=1,iy)
      end do
      do k = kz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((q4(i,j,k),i=1,jx),j=1,iy)
      end do
      noutrec = noutrec + 1
      write (64,rec=noutrec) ((ps4(i,j)+ptop,i=1,jx),j=1,iy)
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
      real :: ptop
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
      write (64,rec=noutrec) idate , jx , iy , kz
      do k = kz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((u4(i,j,k),i=1,jx),j=1,iy)
      end do
      do k = kz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((v4(i,j,k),i=1,jx),j=1,iy)
      end do
      do k = kz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((t4(i,j,k),i=1,jx),j=1,iy)
      end do
      do k = kz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((q4(i,j,k),i=1,jx),j=1,iy)
      end do
      noutrec = noutrec + 1
      write (64,rec=noutrec) ((ps4(i,j)+ptop,i=1,jx),j=1,iy)
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
      real :: ptop
      real , dimension(jx,iy,4) :: qs3 , ti3 , ts3
      real , dimension(jx,iy) :: snow
      intent (in) idate , lsmtyp , ptop , snow , ti3 , ts3 , qs3
!
! Local variables
!
      integer :: i , j , k
!
!     THIS ROUTINE WRITES OUT AN INPUT FILE FOR THE RCM
!
      noutrec = noutrec + 1
      write (64,rec=noutrec) idate , jx , iy , kz
      do k = kz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((u4(i,j,k),i=1,jx),j=1,iy)
      end do
      do k = kz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((v4(i,j,k),i=1,jx),j=1,iy)
      end do
      do k = kz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((t4(i,j,k),i=1,jx),j=1,iy)
      end do
      do k = kz , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((q4(i,j,k),i=1,jx),j=1,iy)
      end do
      noutrec = noutrec + 1
      write (64,rec=noutrec) ((ps4(i,j)+ptop,i=1,jx),j=1,iy)
      noutrec = noutrec + 1
      write (64,rec=noutrec) ts4
 
      if ( lsmtyp=='USGS' ) then
        do k = 1 , 4
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((qs3(i,j,k),i=1,jx),j=1,iy)
        end do
        do k = 1 , 4
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((ti3(i,j,k),i=1,jx),j=1,iy)
        end do
        do k = 1 , 4
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((ts3(i,j,k),i=1,jx),j=1,iy)
        end do
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((snow(i,j),i=1,jx),j=1,iy)
      end if
!
      end subroutine writefs

      end module mod_write
