      program readll
      implicit none
      include 'ERAHI.param'
C
C A      SET PARAMETERS
C
C  X X X X X X X                                      X X X X X X X X X
C  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
C
C  X X X X X   SET 1 : PARAMETERS FOR ERA40 REANALYSIS DATASET   X X X X
C A1
      INTEGER NLATS,NLONS
      PARAMETER (NLATS=160,NLONS=320)
C
C  NLATS  = NUMBER OF LATITUDES  ON ERA40 original GRID.
C  NLONS  = NUMBER OF LONGITUDES ON ERA40 original GRID.
C
      real*4  ps
      common /var_ps/ ps(NLONS,NLATS)
C
      integer bytes_per_kb
      integer bytes_per_mb
      integer bytes_per_wd
      integer words_per_mb 
      integer n80
      integer stp
      
      parameter( bytes_per_kb = 1024 )
      parameter( bytes_per_mb = bytes_per_kb * bytes_per_kb )
      parameter( bytes_per_wd = 4 )
      parameter( words_per_mb = bytes_per_mb / bytes_per_wd )
      parameter( n80 = 80 )
      parameter( stp = 159 )
      real*4 lats
      real*4 area
      common /dimson/ lats(n80),area(4)
      INTEGER IDATE
      INTEGER MDATE
      COMMON /DATENUM/MDATE(67204)
      INTEGER NSTART,NNNEND,NNN
      character*40 finame
      character*14 foname
      character*2 chmon
      INTEGER MYEAR,MONTH
      integer i,j,k
c
      CALL INITDATE
      CALL FINDDATE(NSTART,IDATE1)
      CALL FINDDATE(NNNEND,IDATE2)
      DO NNN=NSTART,NNNEND
         IDATE=MDATE(NNN)
         MYEAR=IDATE/1000000
         MONTH=(IDATE-MYEAR*1000000)/10000
         if(MONTH.eq.1) then
           chmon='01'
         else if(MONTH.eq.2) then
           chmon='02'
         else if(MONTH.eq.3) then
           chmon='03'
         else if(MONTH.eq.4) then
           chmon='04'
         else if(MONTH.eq.5) then
           chmon='05'
         else if(MONTH.eq.6) then
           chmon='06'
         else if(MONTH.eq.7) then
           chmon='07'
         else if(MONTH.eq.8) then
           chmon='08'
         else if(MONTH.eq.9) then
           chmon='09'
         else if(MONTH.eq.10) then
           chmon='10'
         else if(MONTH.eq.11) then
           chmon='11'
         else if(MONTH.eq.12) then
           chmon='12'
         else
           write(*,*) 'MONTH outside of [1,12]'
           stop
         endif
         write(finame,10)MYEAR,chmon,IDATE
 10      format('../DATA/ERAHI/',I4,'/',A2,'/ps.',I10,'.grib')
         write(foname,20) IDATE
 20      format('EHI_',I10)
         open(39,file=foname,form='unformatted',recl=NLONS*NLATS*IBYTE
     &          ,access='direct')
         call get_ps(finame)
         write(39,rec=3) ps
         close(39)
         write(*,*) 'Surface pressure is added for DATE: ',IDATE
      ENDDO
      stop
      end
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE INITDATE
      IMPLICIT NONE
C
      INTEGER MDATE
      COMMON /DATENUM/MDATE(67204)
      INTEGER MBASE,NBASE,NREC,NYEAR,MON,NDAY,I,M
C
      NREC=0
      DO NYEAR=1957,2002
      MBASE = NYEAR*1000000
      DO MON=1,12
         MBASE = MBASE+10000
         IF(MON.EQ.1.OR.MON.EQ.3.OR.MON.EQ.5.OR.MON.EQ.7
     &              .OR.MON.EQ.8.OR.MON.EQ.10.OR.MON.EQ.12) THEN
            NDAY=31
         ELSE IF(MON.EQ.4.OR.MON.EQ.6.OR.MON.EQ.9.OR.MON.EQ.11)THEN
            NDAY=30
         ELSE
            IF(MOD(NYEAR,4).EQ.0) THEN
               NDAY=29
            ELSE
               NDAY=28
            ENDIF
         ENDIF
         NBASE = MBASE
         DO I=1,NDAY
            NBASE = NBASE+100
            DO M=1,4
               NREC=NREC+1
               IF(M.EQ.1) THEN
                  MDATE(NREC)=NBASE
               ELSE IF(M.EQ.2) THEN
                  MDATE(NREC)=NBASE+6
               ELSE IF(M.EQ.3) THEN
                  MDATE(NREC)=NBASE+12
               ELSE
                  MDATE(NREC)=NBASE+18
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      ENDDO
c     WRITE(*,*) 'NREC = ',NREC
      RETURN
      END
C
C  X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
C
      SUBROUTINE FINDDATE(NPOS,IDATE)
      IMPLICIT NONE
      INTEGER MDATE
      COMMON /DATENUM/MDATE(67204)
C
      INTEGER NPOS,IDATE
      INTEGER I
C
      I=0
 10   CONTINUE
      I=I+1
      IF(MDATE(I).EQ.IDATE) THEN
         NPOS=I
         GO TO 200
      ENDIF
      IF(I.GT.67204) GOTO 100
      GO TO 10
 100  WRITE(*,*) 'ERROR IN FINDDATE'
      STOP
 200  RETURN
      END
c23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine get_ps(finame)
c     ------------------------------------------------------------------
      implicit none
      character*40 finame
C
C A      SET PARAMETERS
C
C  X X X X X X X                                      X X X X X X X X X
C  X X X X X X X    USER DEFINED PARAMETERS FOLLOW    X X X X X X X X X
C
C  X X X X X   SET 1 : PARAMETERS FOR ERA40 REANALYSIS DATASET   X X X X
C A1
      INTEGER NLATS,NLONS
      PARAMETER (NLATS=160,NLONS=320)
C
C  NLATS  = NUMBER OF LATITUDES  ON ERA40 original GRID.
C  NLONS  = NUMBER OF LONGITUDES ON ERA40 original GRID.
C
      real*4  ps
      common /var_ps/ ps(NLONS,NLATS)
C
      integer bytes_per_kb
      integer bytes_per_mb
      integer bytes_per_wd
      integer words_per_mb 
      integer n80
      integer stp
      
      parameter( bytes_per_kb = 1024 )
      parameter( bytes_per_mb = bytes_per_kb * bytes_per_kb )
      parameter( bytes_per_wd = 4 )
      parameter( words_per_mb = bytes_per_mb / bytes_per_wd )
      parameter( n80 = 80 )
      parameter( stp = 159 )
      real*4 lats
      real*4 area
      common /dimson/ lats(n80),area(4)

      integer iarraysphh(words_per_mb)

      integer iarrayps(words_per_mb)
      real*8  dave

      real rarraysphh(words_per_mb)
      real rarrayrn80(words_per_mb)

      integer unitsphh

      integer ierror

      integer bytes_read_into_sphh
      integer bytes_read_into_vod

      integer words_read_into_rn80
      integer words_read_into_uv

      integer intv(n80)
      real    realv(n80)
      character*30 charv(n80)

      integer i,j,k

      integer intout, intf
      external intout, intf

c     ------------------------------------------------------------------
c     Open input GRIB file containing spherical harmonic records:

      call pbopen( unitsphh, finame, 'r', ierror )
      call check_error_status('pbopen', ierror)
c     ------------------------------------------------------------------
c     Loop over records:

c     ------------------------------------------------------------------
c     Read 160x320 reduced N80 Gaussian grid records from input GRIB
c     file:

      call pbgrib( unitsphh, iarraysphh, bytes_per_mb,
     &             bytes_read_into_sphh, ierror )

      if ( ierror .eq. -1 ) then
c       end of file (EOF)
        goto 200
      else
        call check_error_status('pbgrib',ierror)
      end if

c     ------------------------------------------------------------------
c     Set up output grid, it being a 160x320 regular N80 Gaussian grid:

      ierror = intout('area',                  intv, area,  charv    )
      call check_error_status('intout', ierror)
      ierror = intout('user_regular_gaussian', n80,   realv, charv   )
      call check_error_status('intout', ierror)
      ierror = intout('g_lats',                intv, lats,  charv    )
      call check_error_status('intout', ierror)
      ierror = intout('truncation',            stp,  realv, charv    )
      call check_error_status('intout', ierror)

c     ------------------------------------------------------------------
c     Compute 160x320 regular N80 Gaussian scalar field from 160x320
c     reduced N80 Gaussian grid:

      words_read_into_rn80 = words_per_mb

      ierror = intf( iarraysphh,  words_per_mb,         rarraysphh,
     &               iarrayps,  words_read_into_rn80, rarrayrn80  )
      call check_error_status('intf', ierror)
      call get_field_data(words_per_mb,iarrayps,rarrayrn80)
      k=0
      do j=NLATS,1,-1
      do i=1,NLONS
         k=k+1
         ps(i,j)=exp(rarrayrn80(k))*0.01
      enddo
      enddo

      dave=0.0d0
      do i=1,NLONS
         dave=dave+ps(i,2)
      enddo
      dave=dave/FLOAT(NLONS)
      do i=1,NLONS
         ps(i,1)=(dave+ps(i,2))*0.5
      enddo

c     ------------------------------------------------------------------
c     Continue statement for loop:

c     ------------------------------------------------------------------
c     Continue after EOF:

200   continue

c     ------------------------------------------------------------------
c     Close input and output GRIB files:

      call pbclose( unitsphh, ierror )
      call check_error_status('pbclose', ierror)

c     ------------------------------------------------------------------
      return
      end
      block data
      implicit none
      integer bytes_per_kb
      integer bytes_per_mb
      integer bytes_per_wd
      integer words_per_mb 
      integer n80
      integer stp
      
      parameter( bytes_per_kb = 1024 )
      parameter( bytes_per_mb = bytes_per_kb * bytes_per_kb )
      parameter( bytes_per_wd = 4 )
      parameter( words_per_mb = bytes_per_mb / bytes_per_wd )
      parameter( n80 = 80 )
      parameter( stp = 159 )
      real*4 lats
      real*4 area
      common /dimson/ lats(n80),area(4)
c     ------------------------------------------------------------------
c     Initialize lats and area, to be used in 'intout':

      data lats/
     + 89.142, 88.029, 86.911, 85.791, 84.670, 83.549, 82.428, 81.307,
     + 80.185, 79.064, 77.943, 76.821, 75.700, 74.578, 73.457, 72.336,
     + 71.214, 70.093, 68.971, 67.850, 66.728, 65.607, 64.485, 63.364,
     + 62.242, 61.121, 60.000, 58.878, 57.757, 56.635, 55.514, 54.392,
     + 53.271, 52.149, 51.028, 49.906, 48.785, 47.663, 46.542, 45.420,
     + 44.299, 43.177, 42.056, 40.934, 39.813, 38.691, 37.570, 36.448,
     + 35.327, 34.205, 33.084, 31.962, 30.841, 29.719, 28.598, 27.476,
     + 26.355, 25.234, 24.112, 22.991, 21.869, 20.748, 19.626, 18.505,
     + 17.383, 16.262, 15.140, 14.019, 12.897, 11.776, 10.654, 9.533,
     + 8.411,  7.290,  6.168,  5.047,  3.925,  2.804,  1.682,  0.561 /

      data area/ 4*0.0 /
      end
c23456789012345678901234567890123456789012345678901234567890123456789012
c     Fortran 77 file/subroutine subr_check_error_status.f
c     Subroutine for checking the error status of integer error codes
c     returned by EMOSLIB routines. A nonzero error code indicates a
c     fatal error: 'check_error_status' will print the ECMWF error
c     message associated with the error code, and will then stop
c     execution of the program. The possible EMOSLIB routines are
c     'pbopen', 'pbwrite', 'pbclose', 'pbgrib', 'intf', 'intuvp',
c     'intin', and 'intout', all character strings. If a subroutine not
c     listed here is used as the first argument, the execution of
c     the program is stopped in this case too.

      subroutine check_error_status( emos_routine, error_code )

      character*(*) emos_routine
      integer error_code

      if ( emos_routine .eq. 'pbopen' ) then
         if ( error_code .ne. 0 ) then
           if (      error_code .eq. -1 ) then
             write(*,'(A80)') 'pbopen: could not open the file'
             stop
           else if ( error_code .eq. -2 ) then
             write(*,'(A80)') 'pbopen: invalid file name'
             stop
           else if ( error_code .eq. -3 ) then
             write(*,'(A80)') 'pbopen: invalid open mode specified'
             stop
           else
             write(*,'(A80)') 'pbopen: unspecified error code' 
             stop
           end if
         else
            return
         end if 
      else if ( emos_routine .eq. 'pbclose' ) then
         if ( error_code .ne. 0 ) then
           if (      error_code .eq. -1 ) then
             write(*,'(A80)') 'pbclose: error closing the file'
             stop
           else
             write(*,'(A80)') 'pbclose: unspecified error code'
             stop
           end if
         else
            return
         end if
      else if ( emos_routine .eq. 'pbwrite' ) then
         if ( error_code .lt. 0 ) then
           if (      error_code .eq. -1 ) then
             write(*,'(A80)') 'pbwrite: error writing to the file'
             stop
           else
             write(*,'(A80)') 'pbwrite: unspecified error code'
             stop
           end if
         else
            return
         end if
      else if ( emos_routine .eq. 'pbgrib' ) then
         if ( error_code .ne. 0 ) then
           if ( error_code .eq. -2 ) then
             write(*,'(A80)') 'pbgrib: error in file-handling'
             write(*,'(A80)') 'pbgrib: e.g. file may contain a trun-'
             write(*,'(A80)') 'cated product                        '
             stop
           else if ( error_code .eq. -3 ) then
             write(*,'(A80)') 'pbgrib: karray is not large enough'
             stop
           else
             write(*,'(A80)') 'pbgrib: unspecified error code'
             stop
           end if
         else
            return
         end if
      else if ( emos_routine .eq. 'intf' ) then
         if ( error_code .ne. 0 ) then
             write(*,'(A80)') 'intf: call not successful'
             stop
         else
            return
         end if
      else if ( emos_routine .eq. 'intuvp' ) then
         if ( error_code .ne. 0 ) then
             write(*,'(A80)') 'intuvp: call not successful'
             stop
         else
            return
         end if
      else if ( emos_routine .eq. 'intout' ) then
         if ( error_code .ne. 0 ) then
             write(*,'(A80)') 'intout: call not successful'
             stop
         else
            return
         end if
      else if ( emos_routine .eq. 'intin' ) then
         if ( error_code .ne. 0 ) then
             write(*,'(A80)') 'intin: call not successful'
             stop
         else
            return
         end if
      else
             write(*,'(A80)') 'check_error_status: ' // emos_routine
             write(*,'(A80)') 'check_error_status: unknown emos routine'
             write(*,'(A80)') 'check_error_status handles pbopen,      '
             write(*,'(A80)') 'pbclose, pbwrite, pbgrib, intf, intuvp, '
             write(*,'(A80)') 'intout, and intin                       '
             stop
      end if

      return
      end
c23456789012345678901234567890123456789012345678901234567890123456789012
c     Fortran 77 file/subroutine: subr_get_field_data.f
c     Assumes an integer GRIB record as input. Uses "gribex" to decode
c     the GRIB record and returns the field data values.
c     ------------------------------------------------------------------
c     IN:
c          words_ioarray: integer, length of integer array containing
c                         GRIB record, in words, and also length of
c                         real output array containing field data values
c                         from decoded input GRIB record
c                 iarray: integer, array containing GRIB record to be
c                         decoded
c     OUT:
c                 rarray: real, array containing field data values
c                         from decoded GRIB record 
c     ------------------------------------------------------------------
      subroutine get_field_data( words_ioarray, iarray, rarray )
c     ------------------------------------------------------------------
c     'Implicit none':
      implicit none

      integer words_ioarray
      integer iarray(words_ioarray)
      real rarray(words_ioarray)

c     ------------------------------------------------------------------
c     Declare variables for use in gribex:
c     ------------------------------------
c     gribex( ksec0, ksec1, ksec2, psec2, ksec3, psec3,
c    +        ksec4, psec4, klenp, kgrib, kleng,
c    +        kword, hoper, kret )

c     ksec0(1): number of 'octets' in GRIB record
c     ksec0(2): GRIB edition number
c     -------------------------------------------
      integer ksec0(2)

c     Product definition section:
c     ---------------------------
      integer ksec1(1024)

c     Grid description section:
c     -------------------------
      integer ksec2(1024)
      real psec2(1024)

c     Bitmap section:
c     ---------------
      integer ksec3(2)
      real psec3(2)

c     Binary data section:
c     --------------------
      integer ksec4(1024)

c     Field data values:
c     ------------------
c     real psec4(words_per_mb)
c     integer klenp
c     parameter( klenp = words_per_mb )
c     (will use real array 'rarray' argument for 'psec4',
c     and  integer words_ioarray argument for 'klenp')

c     Grib formatted binary data:
c     ---------------------------
c     integer kgrib(words_per_mb)
c     integer kleng
c     parameter( kleng = words_per_mb )
c     (will use integer array 'iarray' argument for 'kgrib',
c     and integer words_ioarray argument for 'kleng')

c     Return value from encoding,
c     (number of elements of kgrib occupied by coded data):
c     -----------------------------------------------------
      integer kword
c     (Used only for encoding, but have to supply in any case.)

c     Coding options:
c     ---------------
      character*1 hoper

c     Error handling:
c     ---------------
      integer kret
      kret = 0
c     Abort if an error is encountered. (An initial nonzero value
c     will not cause gribex to abort even if an error is encountered.)
c     Upon return, kret < 0 implies handle warning, kret > 0 implies
c     handle error.

c     ------------------------------------------------------------------
c     Begin 'executable' code:
c     ------------------------------------------------------------------
c     Decode GRIB record:

      hoper = 'D'

      call gribex( ksec0, ksec1, ksec2, psec2, ksec3, psec3,
     +             ksec4, rarray, words_ioarray, iarray, words_ioarray,
     +             kword, hoper, kret )

c     ------------------------------------------------------------------
c     If kret upon return is > 0, print error message and stop program.

      if ( kret > 0 ) then
        write(*,'(A80)') 'gribex: error encountered'
        write(*,'(A80)') 'stopping execution in subroutine ' //
     +                   'get_field_data'
        stop
      else
        return
      end if

c     ------------------------------------------------------------------
      return
      end
