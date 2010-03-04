
SUBROUTINE wrf_message( str )
  IMPLICIT NONE
  CHARACTER*(*) str
#if defined( DM_PARALLEL ) && ! defined( STUBMPI)
  write(0,*) str
#endif
  print*, str
END SUBROUTINE wrf_message

