      program resetsst
      use mod_regcm_param , only : jx , iy
      implicit none
!
! Local variables
!
      integer :: n , nday , nmbr , nmo , nyear
      real , dimension(jx,iy) :: sst
!
      nmbr = 8      ! number need be changed for your own SST
!
      open (10,file='SST.RCM',form='unformatted')
      open (20,file='newSST.RCM',form='unformatted')
      do n = 1 , nmbr
        read (10) nday , nmo , nyear , sst
        print * , 'NDAY=' , nday , ' NMO=' , nmo , ' NYEAR=' , nyear
 
!       write your own code for changing SST if necessary
 
        write (20) nday , nmo , nyear , sst
      end do
!
      end program resetsst
