      module mod_write

      implicit none

      integer :: noutrec

      contains

      subroutine writef(u,v,t,q,px,ts,ptop,ni,nj,nk,idate)
      implicit none
!
! Dummy arguments
!
      integer :: idate , ni , nj , nk
      real :: ptop
      real , dimension(ni,nj) :: px , ts
      real , dimension(ni,nj,nk) :: q , t , u , v
      intent (in) idate , ni , nj , nk , ptop , px , q , t , ts , u , v
!
! Local variables
!
      integer :: i , j , k
!
!     THIS ROUTINE WRITES OUT AN INPUT FILE FOR THE RCM
!
!     PRINT *,'WRITING OUTPUT:  IDATE= ',IDATE
      noutrec = noutrec + 1
      write (64,rec=noutrec) idate , ni , nj , nk
      do k = nk , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((u(i,j,k),i=1,ni),j=1,nj)
      end do
      do k = nk , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((v(i,j,k),i=1,ni),j=1,nj)
      end do
      do k = nk , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((t(i,j,k),i=1,ni),j=1,nj)
      end do
      do k = nk , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((q(i,j,k),i=1,ni),j=1,nj)
      end do
      noutrec = noutrec + 1
      write (64,rec=noutrec) ((px(i,j)+ptop,i=1,ni),j=1,nj)
      noutrec = noutrec + 1
      write (64,rec=noutrec) ts
!
      end subroutine writef

      subroutine writef2(u,v,t,q,px,ts,ptop,ni,nj,nk,idate)
      implicit none
!
! Dummy arguments
!
      integer :: idate , ni , nj , nk
      real :: ptop
      real , dimension(ni,nj,nk) :: q , t , u , v
      real , dimension(ni,nj) :: px , ts
      intent (in) idate , ni , nj , nk , ptop , px , q , t , ts , u , v
!
! Local variables
!
      integer :: i , j , k
!
!     THIS ROUTINE WRITES OUT AN INPUT FILE FOR THE RCM
!
!     PRINT *,'WRITING OUTPUT:  IDATE= ',IDATE
      noutrec = noutrec + 1
      write (64,rec=noutrec) idate , ni , nj , nk
      do k = nk , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((u(i,j,k),i=1,ni),j=1,nj)
      end do
      do k = nk , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((v(i,j,k),i=1,ni),j=1,nj)
      end do
      do k = nk , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((t(i,j,k),i=1,ni),j=1,nj)
      end do
      do k = nk , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((q(i,j,k),i=1,ni),j=1,nj)
      end do
      noutrec = noutrec + 1
      write (64,rec=noutrec) ((px(i,j)+ptop,i=1,ni),j=1,nj)
      noutrec = noutrec + 1
      write (64,rec=noutrec) ts
!
      end subroutine writef2

      subroutine writefs(u,v,t,q,px,ts,qs3,ti3,ts3,snow,ptop,ni,nj,nk,  &
                       & idate,lsmtyp)
      implicit none
!
! Dummy arguments
!
      integer :: idate , ni , nj , nk
      character(4) :: lsmtyp
      real :: ptop
      real , dimension(ni,nj) :: px , snow , ts
      real , dimension(ni,nj,nk) :: q , t , u , v
      real , dimension(ni,nj,4) :: qs3 , ti3 , ts3
      intent (in) idate , lsmtyp , ni , nj , nk , ptop , px , q , qs3 , &
                & snow , t , ti3 , ts , ts3 , u , v
!
! Local variables
!
      integer :: i , j , k
!
!     THIS ROUTINE WRITES OUT AN INPUT FILE FOR THE RCM
!
      noutrec = noutrec + 1
      write (64,rec=noutrec) idate , ni , nj , nk
      do k = nk , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((u(i,j,k),i=1,ni),j=1,nj)
      end do
      do k = nk , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((v(i,j,k),i=1,ni),j=1,nj)
      end do
      do k = nk , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((t(i,j,k),i=1,ni),j=1,nj)
      end do
      do k = nk , 1 , -1
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((q(i,j,k),i=1,ni),j=1,nj)
      end do
      noutrec = noutrec + 1
      write (64,rec=noutrec) ((px(i,j)+ptop,i=1,ni),j=1,nj)
      noutrec = noutrec + 1
      write (64,rec=noutrec) ts
 
      if ( lsmtyp=='USGS' ) then
        do k = 1 , 4
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((qs3(i,j,k),i=1,ni),j=1,nj)
        end do
        do k = 1 , 4
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((ti3(i,j,k),i=1,ni),j=1,nj)
        end do
        do k = 1 , 4
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((ts3(i,j,k),i=1,ni),j=1,nj)
        end do
        noutrec = noutrec + 1
        write (64,rec=noutrec) ((snow(i,j),i=1,ni),j=1,nj)
      end if
!
      end subroutine writefs

      end module mod_write
