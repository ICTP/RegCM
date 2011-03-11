      program tograds2
      implicit none
!
! Local variables
!
      real(4) , dimension(192,96) :: a
      integer :: i , j , k
      integer(2) , dimension(192,96) :: ia
      real(8) :: offset , xscale
!
      open (10,file='/home/RAID2-D13/EH5OM/RF/1961/EHgRF1961JAN',       &
           &form='unformatted',recl=192*96*2+16,access='direct')
      open (20,file='JAN61g.dat',form='unformatted',recl=192*96*4,      &
          & access='direct')
      do k = 1 , 105
        read (10,rec=k) offset , xscale , ia
        do j = 1 , 96
          do i = 1 , 192
            a(i,j) = real(dble(ia(i,97-j))*xscale + offset)
          end do
        end do
        write (20,rec=k) a
      end do
      end program tograds2
