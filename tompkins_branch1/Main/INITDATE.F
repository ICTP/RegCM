      SUBROUTINE INITDATE
      IMPLICIT NONE
C
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
#ifdef MPP1
      include 'mpiregcm.h'
#endif
      INTEGER MDATEZ
      COMMON /DATENUM/MDATEZ(289276)
      INTEGER MBASE,NBASE,NREC,NYEAR,MON,NDAY,I,M
C
      IF(IBDYFRQ.EQ.6) THEN
        NREC=0
        DO NYEAR=1948,2145
          MBASE = NYEAR*1000000
          DO MON=1,12
            MBASE = MBASE+10000
            IF(MON.EQ.1.OR.MON.EQ.3.OR.MON.EQ.5.OR.MON.EQ.7
     &                 .OR.MON.EQ.8.OR.MON.EQ.10.OR.MON.EQ.12) THEN
              NDAY=31
            ELSE IF(MON.EQ.4.OR.MON.EQ.6.OR.MON.EQ.9.OR.MON.EQ.11)THEN
              NDAY=30
            ELSE
              NDAY=28
              IF(MOD(NYEAR,4)  .EQ.0) NDAY=NDAY+1
              IF(MOD(NYEAR,100).EQ.0) NDAY=NDAY-1
              IF(MOD(NYEAR,400).EQ.0) NDAY=NDAY+1
            ENDIF
            NBASE = MBASE
            DO I=1,NDAY
              NBASE = NBASE+100
              DO M=1,4
                NREC=NREC+1
                IF(M.EQ.1) THEN
                  MDATEZ(NREC)=NBASE
                ELSE IF(M.EQ.2) THEN
                  MDATEZ(NREC)=NBASE+6
                ELSE IF(M.EQ.3) THEN
                  MDATEZ(NREC)=NBASE+12
                ELSE
                  MDATEZ(NREC)=NBASE+18
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE IF(IBDYFRQ.EQ.3) THEN
        NREC=0
        DO NYEAR=1948,2046
          MBASE = NYEAR*1000000
          DO MON=1,12
            MBASE = MBASE+10000
            IF(MON.EQ.1.OR.MON.EQ.3.OR.MON.EQ.5.OR.MON.EQ.7
     &                 .OR.MON.EQ.8.OR.MON.EQ.10.OR.MON.EQ.12) THEN
              NDAY=31
            ELSE IF(MON.EQ.4.OR.MON.EQ.6.OR.MON.EQ.9.OR.MON.EQ.11)THEN
              NDAY=30
            ELSE
              NDAY=28
              IF(MOD(NYEAR,4)  .EQ.0) NDAY=NDAY+1
              IF(MOD(NYEAR,100).EQ.0) NDAY=NDAY-1
              IF(MOD(NYEAR,400).EQ.0) NDAY=NDAY+1
            ENDIF
            NBASE = MBASE
            DO I=1,NDAY
              NBASE = NBASE+100
              DO M=1,8
                NREC=NREC+1
                IF(M.EQ.1) THEN
                  MDATEZ(NREC)=NBASE
                ELSE IF(M.EQ.2) THEN
                  MDATEZ(NREC)=NBASE+3
                ELSE IF(M.EQ.3) THEN
                  MDATEZ(NREC)=NBASE+6
                ELSE IF(M.EQ.4) THEN
                  MDATEZ(NREC)=NBASE+9
                ELSE IF(M.EQ.5) THEN
                  MDATEZ(NREC)=NBASE+12
                ELSE IF(M.EQ.6) THEN
                  MDATEZ(NREC)=NBASE+15
                ELSE IF(M.EQ.7) THEN
                  MDATEZ(NREC)=NBASE+18
                ELSE
                  MDATEZ(NREC)=NBASE+21
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE IF(IBDYFRQ.EQ.12) THEN
        NREC=0
        DO NYEAR=1948,2145
          MBASE = NYEAR*1000000
          DO MON=1,12
            MBASE = MBASE+10000
            IF(MON.EQ.1.OR.MON.EQ.3.OR.MON.EQ.5.OR.MON.EQ.7
     &                 .OR.MON.EQ.8.OR.MON.EQ.10.OR.MON.EQ.12) THEN
              NDAY=31
            ELSE IF(MON.EQ.4.OR.MON.EQ.6.OR.MON.EQ.9.OR.MON.EQ.11)THEN
              NDAY=30
            ELSE
              NDAY=28
              IF(MOD(NYEAR,4)  .EQ.0) NDAY=NDAY+1
              IF(MOD(NYEAR,100).EQ.0) NDAY=NDAY-1
              IF(MOD(NYEAR,400).EQ.0) NDAY=NDAY+1
            ENDIF
            NBASE = MBASE
            DO I=1,NDAY
              NBASE = NBASE+100
              DO M=1,2
                NREC=NREC+1
                IF(M.EQ.1) THEN
                  MDATEZ(NREC)=NBASE
                ELSE IF(M.EQ.2) THEN
                  MDATEZ(NREC)=NBASE+12
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        WRITE(*,*) 'Please double check your IBDYFRQ'
        WRITE(*,*) 'We only support IBDYFRQ = 3, 6, or 12'
        WRITE(*,*) 'and we use IBDYFRQ = 6 as the default'
#ifdef MPP1
        WRITE(*,*) 'STOP in INITDATE.f'
        endif
        call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
        STOP
#endif
      ENDIF
#ifdef MPP1
      if (myid .eq. 0) WRITE(*,*) 'NREC = ',NREC
#else
      WRITE(*,*) 'NREC = ',NREC
#endif
      RETURN
      END
