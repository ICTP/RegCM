C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE FINDDATE(NPOS,IDATE)
      IMPLICIT NONE
#ifdef MPP1
      include 'mpiregcm.h'
      integer ierr
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
