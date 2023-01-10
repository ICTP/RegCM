!****************************************************************
!*                                                              *
!*   module MO__SALSA_UPDATE                                 *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to calculate aerosol dynamics                              *
!*                                                              *
!****************************************************************
!----------------------------------------------------------------------------
!!$   Copyright 2014 Atmospheric Research Centre of Eastern Finland,
!!$         Finnish Meteorological Institute, Kuopio, Finland
!!$
!!$   Licensed under the Apache License, Version 2.0 (the "License");
!!$   you may not use this file except in compliance with the License.
!!$   You may obtain a copy of the License at
!!$
!!$       http://www.apache.org/licenses/LICENSE-2.0
!!$
!!$   Unless required by applicable law or agreed to in writing, software
!!$   distributed under the License is distributed on an "AS IS" BASIS,
!!$   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!$   See the License for the specific language governing permissions and
!!$   limitations under the License.
!------------------------------------------------------------------------
!
! -- Added update for cloud bins (05/2014 J Tonttila, FMI)
!
!****************************************************************

MODULE mo_salsa_update

CONTAINS

  SUBROUTINE distr_update(kproma, kbdim, klev,   &
                          paero,  pcloud, pprecp )

    USE mo_submctl
    USE mo_kind, ONLY : dp

    IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid points 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical levels

    TYPE(t_section), INTENT(inout) :: pcloud(kbdim,klev,ncld), & ! Cloud size distribution and properties
                                      pprecp(kbdim,klev,nprc), & ! Rain drop size distribution and properties
                                      paero(kbdim,klev,fn2b)

    !-- Local variables ----------------------
    INTEGER :: ii, jj, kk, mm
    REAL(dp) :: zvpart, znfrac, zvfrac, zVrat, zVilo, zVihi, zVexc, zvdec
    LOGICAL  :: within_bins
    INTEGER :: count
    
    DO jj = 1,klev
       DO ii = 1,kproma
          
          ! ------------------------------------------------------------------------
          ! ************* AEROSOLS **************
          ! ------------------------------------------------------------------------

          within_bins = .FALSE.
          !-- Check if the volume of the bin is within bin limits after update
          count = 0
          DO WHILE(.NOT.within_bins)
             within_bins = .TRUE.

             DO kk = fn2b-1,in1a,-1
                mm = 0
                IF (paero(ii,jj,kk)%numc > nlim) THEN

                   IF (kk == fn2a) CYCLE 

                   ! Dry volume
                   zvpart = sum(paero(ii,jj,kk)%volc(1:7))/paero(ii,jj,kk)%numc

                   ! Smallest bin cannot decrease
                   IF (paero(ii,jj,kk)%vlolim > zvpart .AND. kk == in1a) CYCLE

                   ! Decreasing bins
                   IF(paero(ii,jj,kk)%vlolim > zvpart) THEN
                      mm = kk - 1
                      IF(kk == in2b) mm = fn1a
                      
                      paero(ii,jj,mm)%numc = paero(ii,jj,mm)%numc + paero(ii,jj,kk)%numc 
                      paero(ii,jj,kk)%numc = 0._dp
                      paero(ii,jj,mm)%volc(:) = paero(ii,jj,mm)%volc(:) + paero(ii,jj,kk)%volc(:)
                      paero(ii,jj,kk)%volc(:) = 0._dp
                      CYCLE
                   END IF

                   !-- If size bin has not grown, cycle
                   IF(zvpart <= pi6*paero(ii,jj,kk)%dmid**3) CYCLE

                   !-- volume ratio of the size bin
                   zVrat = paero(ii,jj,kk)%vhilim/paero(ii,jj,kk)%vlolim
                
                   !-- particle volume at the low end of the bin
                   zVilo = 2._dp*zvpart/(1._dp + zVrat)

                   !-- particle volume at the high end of the bin
                   zVihi = zVrat * zVilo
                   
                   !-- volume in the grown bin which exceeds 
                   !   the bin upper limit
                   zVexc = 0.5_dp*(zVihi + paero(ii,jj,kk)%vhilim)

                   !-- number fraction to be moved to the larger bin
                   znfrac = min(1._dp,(zVihi-paero(ii,jj,kk)%vhilim) / (zVihi - zVilo))
          
                   !-- volume fraction to be moved to the larger bin
                   !zvfrac = znfrac * zVexc / (1._dp/2._dp*(zVihi+zVilo))
                   zvfrac = MIN(0.99_dp,znfrac*zVexc/zvpart)

                   IF(zvfrac > 1._dp) WRITE(*,*) 'VIRHE aerosol'
                   IF(zvfrac < 0._dp) WRITE(*,*) 'VIRHE aerosol 0'
                   !-- update bin
                   mm = kk+1
                   !-- volume
                   paero(ii,jj,mm)%volc(:) = paero(ii,jj,mm)%volc(:) &
                        + znfrac * paero(ii,jj,kk)%numc * zVexc * paero(ii,jj,kk)%volc(:) / &
                        sum(paero(ii,jj,kk)%volc(1:7))

                   paero(ii,jj,kk)%volc(:) = paero(ii,jj,kk)%volc(:) &
                        - znfrac * paero(ii,jj,kk)%numc * zVexc * paero(ii,jj,kk)%volc(:) / &
                        sum(paero(ii,jj,kk)%volc(1:7))

                   !-- number
                   paero(ii,jj,mm)%numc = paero(ii,jj,mm)%numc + znfrac * paero(ii,jj,kk)%numc

                   paero(ii,jj,kk)%numc = paero(ii,jj,kk)%numc * (1._dp - znfrac)

                END IF ! nlim

                IF ( paero(ii,jj,kk)%numc > nlim ) THEN
                   zvpart = sum(paero(ii,jj,kk)%volc(1:7))/paero(ii,jj,kk)%numc
                   IF(zvpart > paero(ii,jj,kk)%vhilim) within_bins = .FALSE.
                   IF(zvpart > paero(ii,jj,kk)%vhilim) THEN
                      write(6,*) kk
                   END IF
                   
                END IF

             END DO ! - kk

             count = count + 1
             IF (count > 100) THEN
                WRITE(*,*) 'WARNING: Aerosol bin update not converged'
                EXIT
             END IF

          END DO ! - within bins

          ! ------------------------------------------------------------------------
          ! ************* CLOUD DROPLETS  **************
          ! ------------------------------------------------------------------------
          
          within_bins = .FALSE.
          count = 0
          DO WHILE (.NOT. within_bins) 
             within_bins = .TRUE.

             DO kk = ncld,ica%cur,-1
                mm = 0

                IF ( pcloud(ii,jj,kk)%numc > nlim .AND. sum(pcloud(ii,jj,kk)%volc(1:5)) > 1.e-30_dp ) THEN
                
                   ! Don't convert cloud or rain droplets to anything else here.
                   zvpart = sum(pcloud(ii,jj,kk)%volc(1:5))/pcloud(ii,jj,kk)%numc
                
                   !-- volume ratio of the size bin
                   zVrat = pcloud(ii,jj,kk)%vhilim/pcloud(ii,jj,kk)%vlolim
                
                   !-- particle volume at the low end of the bin
                   zVilo = 2._dp*zvpart/(1._dp + zVrat)

                   !-- particle volume at the high end of the bin
                   zVihi = zVrat * zVilo

                   !-- Decreasing droplets
                   IF ( zvpart < pi6*pcloud(ii,jj,kk)%vlolim .AND.  &
                        (kk /= ica%cur .AND. kk /= icb%cur)    ) THEN 

                      !-- Volume in the decreased bin which is below the bin lower limit
                      zVexc = 0.5_dp*(zVilo + pcloud(ii,jj,kk)%vlolim)
                      
                      !-- Number fraction to be moved to the smaller bin
                      znfrac = min(1._dp,(pcloud(ii,jj,kk)%vlolim-zVilo) / (zVihi-zVilo))
                      
                      !-- Index for the smaller bin
                      mm = kk - 1

                   !-- Increasing droplets
                   ELSE IF ( zvpart > pi6*pcloud(ii,jj,kk)%dmid**3 .AND.  &
                             (kk /= fca%cur .AND. kk /= icb%cur)     )  THEN  ! Increasing droplets  

                      !-- volume in the grown bin which exceeds the bin upper limit
                      zVexc = 0.5_dp*(zVihi + pcloud(ii,jj,kk)%vhilim)

                      !-- number fraction to be moved to the larger bin
                      znfrac = min(1._dp,(zVihi-pcloud(ii,jj,kk)%vhilim) / (zVihi-zVilo))
          
                      !-- Index for the larger bin
                      mm = kk + 1

                      IF (mm > 14) WRITE(*,*) 'update, heihoi'

                   ELSE  ! Particle size unchanged
                      CYCLE
                   END IF

                   IF ( mm > 14) WRITE(*,*) 'update', kk, fca%cur

                   !-- volume fraction to be moved 
                   !zvfrac = znfrac * zVexc / (0.5_dp*(zVihi+zVilo))
                   !zvfrac = min(zvfrac,1._dp)
                   zvfrac = MIN(0.99_dp,znfrac*zVexc/zvpart)
                   IF(zvfrac > 1._dp) THEN
                      WRITE(*,*) 'VIRHE cloud'
                      WRITE(*,*) zvfrac,zvpart,pi6*pcloud(ii,jj,kk)%dmid**3
                   END IF

                   IF(zvfrac < 0._dp) WRITE(*,*) 'VIRHE cloud 0'
                   !-- volume
                   pcloud(ii,jj,mm)%volc(:) = pcloud(ii,jj,mm)%volc(:)     &
                        + zvfrac*pcloud(ii,jj,kk)%volc(:)
                   
                   pcloud(ii,jj,kk)%volc(:) = pcloud(ii,jj,kk)%volc(:)     &
                        - zvfrac*pcloud(ii,jj,kk)%volc(:)
 
                   !-- number
                   pcloud(ii,jj,mm)%numc = pcloud(ii,jj,mm)%numc    &
                        + znfrac*pcloud(ii,jj,kk)%numc
                   
                   pcloud(ii,jj,kk)%numc = pcloud(ii,jj,kk)%numc    &
                        - znfrac*pcloud(ii,jj,kk)%numc
              
                END IF !nlim
                
                IF ( pcloud(ii,jj,kk)%numc > nlim .AND.  sum(pcloud(ii,jj,kk)%volc(1:5)) > 1.e-30_dp ) THEN
                   zvpart = sum(pcloud(ii,jj,kk)%volc(1:5))/pcloud(ii,jj,kk)%numc

                   IF(zvpart > pcloud(ii,jj,kk)%vhilim  .OR.   &
                      zvpart < pcloud(ii,jj,kk)%vlolim) THEN
                      within_bins = .FALSE.
                   END IF

                END IF

             END DO !kk

             count = count + 1
             IF (count > 100) THEN
                WRITE(*,*) 'WARNING: Cloud bin update not converged'
                WRITE(*,*) pcloud(ii,jj,1:7)%numc
                WRITE(*,*) pcloud(ii,jj,1:7)%volc(1)
                WRITE(*,*) pcloud(ii,jj,1:7)%volc(2)
                WRITE(*,*) pcloud(ii,jj,1:7)%volc(3)
                WRITE(*,*) pcloud(ii,jj,1:7)%volc(4)
                WRITE(*,*) pcloud(ii,jj,1:7)%volc(5)
                WRITE(*,*) pcloud(ii,jj,1:7)%volc(6)
                WRITE(*,*) pcloud(ii,jj,1:7)%volc(7)
                WRITE(*,*) pcloud(ii,jj,1:7)%volc(8)
                EXIT
             END IF

          END DO !within_bins


          ! ------------------------------------------------------------------------
          ! ************* RAIN DROPS **************
          ! Everything else the same as with cloud 
          ! droplets & aerosols, except that the rain 
          ! bins are organized according to the wet 
          ! diameter.
          ! ------------------------------------------------------------------------
          
          within_bins = .FALSE.
          count = 0
          DO WHILE (.NOT. within_bins) 
             within_bins = .TRUE.

             ! -- Juha: now the same for cloud bins
             DO kk = nprc,ira,-1
                mm = 0
                IF ( pprecp(ii,jj,kk)%numc > prlim ) THEN
                
                   zvpart = sum(pprecp(ii,jj,kk)%volc(1:8))/pprecp(ii,jj,kk)%numc
                
                   !-- volume ratio of the size bin
                   zVrat = pprecp(ii,jj,kk)%vhilim/pprecp(ii,jj,kk)%vlolim
                
                   !-- particle volume at the low end of the bin
                   zVilo = 2._dp*zvpart/(1._dp + zVrat)

                   !-- particle volume at the high end of the bin
                   zVihi = zVrat * zVilo

                   ! Calculate the threshold particle volume for decreasing
                   zvdec = (pi6*pprecp(ii,jj,kk)%dmid**3) - &
                           0.2_dp*((pi6*pprecp(ii,jj,kk)%dmid**3) - pprecp(ii,jj,kk)%vlolim)

                   !-- Decreasing droplets - This is now more critical since we are following the wet diameter!!!
                   IF ( zvpart < zvdec .AND. kk /= ira ) THEN 

                      !-- Volume in the decreased bin which is below the bin lower limit
                      zVexc = 0.5_dp*(zVilo + pprecp(ii,jj,kk)%vlolim)
                      
                      !-- Number fraction to be moved to the smaller bin
                      znfrac = min(1._dp,(pprecp(ii,jj,kk)%vlolim-zVilo) / (zVihi-zVilo))
                      IF (znfrac < 0._dp) THEN
                         WRITE(*,*) 'VIRHE, numc precp 0, DEC'
                         write(*,*) zVihi-zVilo
                      END IF

                      !-- Index for the smaller bin
                      mm = kk - 1

                   !-- Increasing droplets
                   ELSE IF ( zvpart > pi6*pprecp(ii,jj,kk)%dmid**3 .AND. kk /= fra )  THEN  ! Increasing droplets  

                      !-- volume in the grown bin which exceeds the bin upper limit
                      zVexc = 0.5_dp*(zVihi + pprecp(ii,jj,kk)%vhilim)

                      !-- number fraction to be moved to the larger bin
                      znfrac = min(.99_dp,(zVihi-pprecp(ii,jj,kk)%vhilim) / (zVihi-zVilo))
          
                      IF (znfrac < 0._dp) THEN
                         WRITE(*,*) 'VIRHE, numc precp 0, INC'
                         write(*,*) zVihi-zVilo
                         WRITE(*,*) zVihi, pprecp(ii,jj,kk)%vhilim,zvpart
                         WRITE(*,*) pi6*pprecp(ii,jj,kk)%dmid**3
                      END IF

                      !-- Index for the larger bin
                      mm = kk + 1

                   ELSE  ! Particle size unchanged
                      CYCLE
                   END IF

                   !-- volume fraction to be moved 
                   !zvfrac = min(0.99,znfrac * zVexc / (0.5_dp*(zVihi+zVilo)))
                   zvfrac = MIN(0.99_dp,znfrac*zVexc/zvpart)

                   IF(zvfrac > 1._dp) WRITE(*,*) 'VIRHE, precp'
                   IF(zvfrac < 0._dp) WRITE(*,*) 'VIRHE, precp 0'

                   !-- volume
                   pprecp(ii,jj,mm)%volc(:) = pprecp(ii,jj,mm)%volc(:)     &
                        + zvfrac*pprecp(ii,jj,kk)%volc(:)
                   
                   pprecp(ii,jj,kk)%volc(:) = pprecp(ii,jj,kk)%volc(:)     &
                        - zvfrac*pprecp(ii,jj,kk)%volc(:)
 
                   !-- number
                   pprecp(ii,jj,mm)%numc = pprecp(ii,jj,mm)%numc    &
                        + znfrac*pprecp(ii,jj,kk)%numc
                   
                   pprecp(ii,jj,kk)%numc = pprecp(ii,jj,kk)%numc    &
                        - znfrac*pprecp(ii,jj,kk)%numc
              
                END IF !nlim
                
                IF ( pprecp(ii,jj,kk)%numc > prlim ) THEN
                   zvpart = sum(pprecp(ii,jj,kk)%volc(:))/pprecp(ii,jj,kk)%numc

                   IF(zvpart > pprecp(ii,jj,kk)%vhilim  .OR.   &
                      zvpart < pprecp(ii,jj,kk)%vlolim) THEN
                      within_bins = .FALSE.
                   END IF

                END IF

             END DO !kk

             count = count + 1
             IF (count > 100) THEN
                WRITE(*,*) 'WARNING: precipitation bin update not converged'
                WRITE(*,*) pprecp(ii,jj,ira:fra)%numc
                WRITE(*,*) pprecp(ii,jj,ira:fra)%volc(8)
                EXIT
             END IF

          END DO !within_bins

       END DO    ! - ii
    END DO       ! - jj


665 FORMAT(999(E11.4,1X))

  END SUBROUTINE distr_update

END MODULE mo_salsa_update
