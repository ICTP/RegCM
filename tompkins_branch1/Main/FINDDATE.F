
      SUBROUTINE FINDDATE(NPOS,IDATE)
      IMPLICIT NONE
#ifdef MPP1
      include 'mpiregcm.h'
#endif
      INTEGER         MDATEZ
      COMMON /DATENUM/MDATEZ(289276)
C
      INTEGER NPOS,IDATE
      INTEGER I
C
      I=0
 10   CONTINUE
      I=I+1
      IF(MDATEZ(I).EQ.IDATE) THEN
         NPOS=I
         GO TO 200
      ENDIF
      IF(I.GT.289276) GOTO 100
      GO TO 10
 100  CONTINUE
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      WRITE(*,*) 'ERROR IN FINDDATE',NPOS,IDATE
#ifdef MPP1
      endif
      call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
      STOP
#endif
 200  RETURN
      END
