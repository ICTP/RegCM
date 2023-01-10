MODULE mo_salsa_sizedist
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
!----------------------------------------------------------------------------
CONTAINS 

  SUBROUTINE size_distribution(kproma, kbdim, klev, &
       n, dpg, sigmag, naero)

    USE mo_submctl, ONLY :      &
         nreg,                      &
         pi6,                       &
         pi,                        &
         in1a,                      &
         fn2a,                      &
         fn2b

    USE mo_salsa_driver, ONLY : aero ! This is needed for size bins spacings

    USE mo_kind, ONLY : dp

    IMPLICIT NONE

    INTEGER, PARAMETER :: nmod = 7

    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid points 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical levels 

    REAL(dp), INTENT(IN) ::         &
         n(nmod)                  , & ! total concentration of a mode
         dpg(nmod)                , & ! geometric-mean diameter of a mode
         sigmag(nmod)                 ! standard deviation of a mode

    REAL(dp), INTENT(OUT) ::        &
         naero(kbdim,klev,fn2b)      ! number concentration  [#/m3]

    !-- local variables
    REAL(dp) ::                     &
         deltadp                      ! bin width [m]

    REAL(dp) :: d1,d2,delta_d,dmid

    INTEGER :: ii, jj, kk,ib

    naero = 0._dp

    DO jj = 1,klev    ! vertical grid
       DO ii = 1,kproma ! horizontal grid

          DO kk = in1a, fn2a
             naero(ii,jj,kk) = 0.0

             d1 = (aero(ii,jj,kk)%vlolim/pi6)**(1._dp/3._dp)
             d2 = (aero(ii,jj,kk)%vhilim/pi6)**(1._dp/3._dp)
             delta_d=(d2-d1)/10
             DO ib = 1,10
                d1=(aero(ii,jj,kk)%vlolim/pi6)**(1._dp/3._dp)+(ib-1.)*delta_d
                d2=d1+delta_d
                dmid=(d1+d2)/2
                deltadp = log(d2/d1)
                
                !deltadp = (aero(ii,jj,kk)%vhilim**(1._dp/3._dp)-aero(ii,jj,kk)%vlolim**(1._dp/3._dp))/   &
                !     pi6**(1._dp/3._dp)
                
                !-- size distribution
                !   ntot = total number, total area, or total volume concentration
                !   dpg = geometric-mean number, area, or volume diameter
                !   n(kk) = number, area, or volume concentration in a bin
!                naero(ii,jj,kk) = naero(ii,jj,kk)+sum(n*deltadp/                        &
!                     (sqrt(2._dp*pi)*log(sigmag))*                 &
!                     exp(-log(aero(ii,jj,kk)%dmid/dpg)**2/(2._dp*log(sigmag)**2)))
                  naero(ii,jj,kk) = naero(ii,jj,kk)+sum(n*deltadp/                        &
                     (sqrt(2._dp*pi)*log(sigmag))*                 &
                     exp(-log(dmid/dpg)**2/(2._dp*log(sigmag)**2)))
                

 !                write(*,*) naero(ii,jj,kk)*1e-6,d1,d2,dmid
             END DO
 !            write(*,*) dpg(1:2),n(1:2)*1e-6
 !            write(*,*) kk,naero(ii,jj,kk)*1e-6,aero(ii,jj,kk)%dmid, (aero(ii,jj,kk)%vlolim/pi6)**(1._dp/3._dp), &
 !                 (aero(ii,jj,kk)%vhilim/pi6)**(1._dp/3._dp)
 !            write(*,*) sum(naero(ii,jj,:))
 !            write(*,*) 'here'
 !            pause
          END DO

       END DO

    END DO

  END SUBROUTINE size_distribution

END MODULE mo_salsa_sizedist
