!****************************************************************
!*	                                                        *
!*   module MO__SALSA_PROPERTIES                                 *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to calculate particle properties during simulation         *
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
!----------------------------------------------------------------------------
MODULE mo_salsa_properties


CONTAINS

  ! fxm: should sea salt form a solid particle when prh is very low
  !  (even though it could be mixed with e.g. sulphate)?
  ! fxm: crashes if no sulphate or sea salt
  ! fxm: do we really need to consider Kelvin effect for regime 2
  !********************************************************************
  !
  ! subroutine WETSIZE()
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculates ambient sizes of particles by equilibrating
  !  soluble fraction of particles with water
  !
  !
  ! Method:
  ! ------- 
  ! Following chemical components are assumed water-soluble
  ! - (ammonium) sulphate (100%)
  ! - sea salt (100 %)
  ! - organic carbon (epsoc * 100%)
  !
  ! Exact thermodynamic considerations neglected
  ! - If particles contain no sea salt, calculation according to
  !  sulphate properties
  ! - If contain sea salt but no sulphate, calculation according to
  !  sea salt properties
  ! - If contain both sulphate and sea salt
  !  -> the molar fraction of these compounds determines
  !     which one of them is used as the basis of calculation
  !     
  ! If sulphate and sea salt coexist in a particle,
  !   it is assumed that the Cl is replaced by sulphate;
  !   thus only either sulphate + organics or sea salt + organics
  !   is included in the calculation of soluble fraction.
  !
  ! Molality parameterizations taken from table 1 of
  !  Tang: Mixed-salt aerosols of atmospheric importance,
  !   JGR, 102 (D2), 1883-1893 (1997)
  !
  !
  ! Interface:
  ! ----------
  ! Called from main aerosol model
  ! 
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005 
  ! Harri Kokkola (FMI) 2006
  ! Matti Niskanen(FMI) 2012
  ! Anton Laakso  (FMI) 2013
  !
  ! Modified for the new aerosol datatype,
  ! Juha Tonttila (FMI) 2014
  !---------------------------------------------------------------------

  SUBROUTINE equilibration(kproma, kbdim, klev,    &
                           prh, ptemp, paero, init )

    USE mo_submctl, ONLY : &
         t_section,    &
         pi6,          & ! pi/6
         in1a, fn1a,   &
         in2a, fn2a,   &
         in2b, fn2b,   &
         boltz,        & ! Boltzmann constant [J/K]
         nlim,         & ! lowest possible particle conc. in a bin [#/m3]
    
                         ! molar masses [kg/mol]
         mwa,          & ! water
                         ! molecular volumes [m3]
         mvsu,         & ! sulphate 
                         ! density [kg/m3]
         rhowa,        & ! water
         
         surfw0,       & ! surface tension of water [J/m2]
         epsoc,        & ! fxm

         rhosu, msu,   & ! properties of compounds
         rhooc, moc,   &
         rhobc, mbc,   &
         rhoss, mss,   &
         rhono, mno,   &
         rhonh, mnh


    USE mo_kind, ONLY : dp

    IMPLICIT NONE

    !-- input variables -------------
    INTEGER, INTENT(in) ::          &
         kproma,                    & ! number of horiz. grid kproma 
         kbdim,                     & ! dimension for arrays 
         klev                         ! number of vertical levels 

    REAL(dp), INTENT(in) ::        &     
         prh(kbdim,klev),          & ! relative humidity [0-1]
         ptemp(kbdim,klev)           ! temperature [K]

    LOGICAL, INTENT(in) :: init  ! TRUE: Initialization call
                                 ! FALSE: Normal runtime: update water content only for 1a

    !-- output variables -------------
    TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,fn2b)     



    !-- local variables --------------
    INTEGER :: ii, jj, kk             ! loop indices
    INTEGER :: count

    REAL(dp) ::      &
         zbinmol(7), &   ! binary molality of individual components [mol/kg]
         zvpart(7),  &   ! volume of chem. compounds in one particle [fxm]
         zke,        &   ! Kelvin term
         zaw,        &   ! water activity [0-1]
         zlwc,       &   ! liquid water content [kg/m3-air]
         zdold,      &   !
         zrh
    REAL(dp) :: zcore,  &
                zdwet


    !----------------------------------------------------------------------
    !-- 1) Regime 1: sulphate and partly water-soluble OC -----------------
    !                This is done for every CALL
    zke = 1.001_dp
    !pdwet = 0._dp
    DO kk = in1a,fn1a      ! size bin
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid

             !-- initialize
             zbinmol = 0._dp
             zdold = 1._dp

             IF ((paero(ii,jj,kk)%numc > nlim)) THEN

                !-- volume of sulphate and OC in one particle [fxm]
                !   + NO/NH

                zvpart(1:2) = paero(ii,jj,kk)%volc(1:2)/paero(ii,jj,kk)%numc
                zvpart(6:7) = paero(ii,jj,kk)%volc(6:7)/paero(ii,jj,kk)%numc ! NO + NH
                
                !-- total volume of one dry particle [fxm] 
                zcore = sum(zvpart(1:2))
                zdwet = paero(ii,jj,kk)%dwet
                
                ! Relative Humidity:
                zrh = prh(ii,jj)
                zrh = MAX(zrh , 0.05_dp)
                zrh = MIN(zrh , 0.98_dp)

                count = 0
                DO WHILE(abs(zdwet/zdold-1.) > 1.e-2_dp) 
                   zdold = max(zdwet,1.e-20_dp)

                   zaw = zrh/zke

                   !-- binary molalities [mol/kg]
                   zbinmol(1) =                      & ! sulphate
                        + 1.1065495e+2_dp            & 
                        - 3.6759197e+2_dp * zaw      &  
                        + 5.0462934e+2_dp * zaw**2   &
                        - 3.1543839e+2_dp * zaw**3   &
                        + 6.770824e+1_dp  * zaw**4 

                   zbinmol(2) = 1./(zaw*mwa)-1./mwa ! organic carbon

                   zbinmol(6) =                      & ! nitric acid
                        + 2.306844303e+1_dp          &
                        - 3.563608869e+1_dp * zaw    &
                        - 6.210577919e+1_dp * zaw**2 &
                        + 5.510176187e+2_dp * zaw**3 &
                        - 1.460055286e+3_dp * zaw**4 &
                        + 1.894467542e+3_dp * zaw**5 &
                        - 1.220611402e+3_dp * zaw**6 &
                        + 3.098597737e+2_dp * zaw**7

                   !
                   ! Calculate the liquid water content (kg/m3-air) using ZSR
                   ! (see e.g. equation (9.98) in Seinfeld and Pandis (1998))
                   !
                   zlwc = (paero(ii,jj,kk)%volc(1)*(rhosu/msu))/zbinmol(1)       + &
                          epsoc * paero(ii,jj,kk)%volc(2)*(rhooc/moc)/zbinmol(2) + &
                          (paero(ii,jj,kk)%volc(6)*(rhono/mno))/zbinmol(6)
                   
                   !-- particle wet radius [m] 
                   zdwet = (zlwc/paero(ii,jj,kk)%numc/rhowa/pi6 + &
                        SUM(zvpart(6:7))/pi6 + zcore/pi6)**(1._dp/3._dp)

                   zke = exp(2._dp*surfw0*mvsu/(boltz*ptemp(ii,jj)*zdwet))
                   !-- Kelvin effect 
                   
                   count = count + 1
                   IF (count > 1000) THEN
                      WRITE(*,*) 'SALSA properties: no convergence!!'
                      EXIT
                   END IF

                END DO

                ! Instead of lwc, use the volume concentration of water from now on 
                ! (easy to convert...)
                !plwc(ii,jj,kk) = zlwc
                paero(ii,jj,kk)%volc(8) = zlwc/rhowa
                
                ! If this is initialization, update the core and wet diameter
                IF (init) THEN
                   paero(ii,jj,kk)%dwet = zdwet
                   paero(ii,jj,kk)%core = zcore
                END IF


             ELSE
                ! If initialization
                !-- 1.2) empty bins given bin average values ----------------- 
                IF (init) THEN
                   paero(ii,jj,kk)%dwet = paero(ii,jj,kk)%dmid
                   paero(ii,jj,kk)%core = pi6*paero(ii,jj,kk)%dmid**3
                END IF
             END IF
          END DO
       END DO
    END DO

    !-- 2) Regime 2a: sulphate, OC, BC and sea salt ----------------------------
    !                 This is done only for initialization call, otherwise the 
    !                 water contents are computed via condensation

    IF (init) THEN

       ! loops over:
       DO kk = in2a,fn2b      ! size bin
          DO jj = 1,klev      ! vertical grid
             DO ii = 1,kproma ! horizontal grid

                zke = 1.02_dp
            
                !-- initialize
                zbinmol = 0._dp
                zdold = 1._dp

                !-- 1) particle properties calculated for non-empty bins ---------
                IF ((paero(ii,jj,kk)%numc > nlim)) THEN

                   !-- volume in one particle [fxm]
                   zvpart = 0._dp
                   zvpart(1:7) = paero(ii,jj,kk)%volc(1:7)/paero(ii,jj,kk)%numc
                   
                   !-- total volume of one dry particle [fxm] 
                   zcore = sum(zvpart(1:5))
                   zdwet = paero(ii,jj,kk)%dwet

                   ! Relative Humidity:
                   zrh = prh(ii,jj)
                   zrh = MAX(zrh , 0.37_dp)
                   zrh = MIN(zrh , 0.98_dp)

                   count = 0
                   DO WHILE(abs(zdwet/zdold-1.) > 1.e-12_dp)
                      zdold = max(zdwet,1.e-20_dp)

                      zaw = zrh/zke

                      !-- binary molalities [mol/kg]
                      zbinmol(1) =                     &  ! sulphate
                           + 1.1065495e+2_dp           & 
                           - 3.6759197e+2_dp * zaw     &  
                           + 5.0462934e+2_dp * zaw**2  &
                           - 3.1543839e+2_dp * zaw**3  &
                           + 6.770824e+1_dp  * zaw**4 
                      
                      zbinmol(2) = 1._dp/(zaw*mwa)-1._dp/mwa ! organic carbon

                      zbinmol(6) =                      & ! nitric acid
                           + 2.306844303e+1_dp          &
                           - 3.563608869e+1_dp * zaw    &
                           - 6.210577919e+1_dp * zaw**2 &
                           + 5.510176187e+2_dp * zaw**3 &
                           - 1.460055286e+3_dp * zaw**4 &
                           + 1.894467542e+3_dp * zaw**5 &
                           - 1.220611402e+3_dp * zaw**6 &
                           + 3.098597737e+2_dp * zaw**7

                      zbinmol(5) =                     &  ! sea salt (NaCl)
                           + 5.875248e+1_dp            &  ! 
                           - 1.8781997e+2_dp * zaw     &  
                           + 2.7211377e+2_dp * zaw**2  &
                           - 1.8458287e+2_dp * zaw**3  &
                           + 4.153689e+1_dp  * zaw**4 
                      
                      !-- calculate the liquid water content (kg/m3-air)
                      zlwc = (paero(ii,jj,kk)%volc(1)*(rhosu/msu))/zbinmol(1) +                 &
                           epsoc * (paero(ii,jj,kk)%volc(2)*(rhooc/moc))/zbinmol(2) +         &
                           (paero(ii,jj,kk)%volc(6)*(rhono/mno))/zbinmol(6)         +         &
                           (paero(ii,jj,kk)%volc(5)*(rhoss/mss))/zbinmol(5)
                      
                      !-- particle wet radius [m] 
                      zdwet = (zlwc/paero(ii,jj,kk)%numc/rhowa/pi6 +  &
                           SUM(zvpart(6:7))/pi6 + zcore/pi6)**(1._dp/3._dp)

                      !-- Kelvin effect 
                      zke = exp(2._dp*surfw0*mvsu/(boltz*ptemp(ii,jj)*zdwet))

                      count = count + 1
                      IF (count > 1000) THEN
                         WRITE(*,*) 'SALSA properties: no convergence!!'
                         EXIT
                      END IF
                      
                   END DO

                   ! Liquid water content; instead of LWC use the volume concentration
                   !plwc(ii,jj,kk)=zlwc
                   paero(ii,jj,kk)%volc(8) = zlwc/rhowa
                   paero(ii,jj,kk)%dwet = zdwet
                   paero(ii,jj,kk)%core = zcore

                ELSE
                   !-- 2.2) empty bins given bin average values ------------------------- 
                   paero(ii,jj,kk)%dwet = paero(ii,jj,kk)%dmid
                   paero(ii,jj,kk)%core = pi6*paero(ii,jj,kk)%dmid**3
                END IF
                
             END DO
          END DO
       END DO

    END IF

  END SUBROUTINE equilibration

END MODULE mo_salsa_properties
