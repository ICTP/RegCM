!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_che_drydep
  !
  ! Chemical and aerosol surface emission and dry deposition
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_che_common
  use mod_che_dust
  use mod_mpmessage
  use mod_service
  use mod_che_ncio
  use mod_che_mppio
  use mod_che_indices

  implicit none

  private

  public :: drydep_aero , drydep_gas , aerodyresis
  public :: a1 , a2 , a3 , c1 , c2 , c3 , c4 , aa1 , aa2 , aa3
  !
  ! Dynamic Viscosity Parameters
  !
  real(rkx) , parameter :: a1 = 1.458e-6_rkx
  real(rkx) , parameter :: a2 = 1.5_rkx
  real(rkx) , parameter :: a3 = 110.4_rkx
  !
  ! Molecular Free Path calculation parameters
  !
  real(rkx) , parameter :: c1 = 6.54e-8_rkx
  real(rkx) , parameter :: c2 = 1.818e-5_rkx
  real(rkx) , parameter :: c3 = 1.013e5_rkx
  real(rkx) , parameter :: c4 = 293.15_rkx
  !
  ! Cunningham slip correction factor parameters
  !
  real(rkx) , parameter :: aa1 = 1.257_rkx
  real(rkx) , parameter :: aa2 = 0.4_rkx
  real(rkx) , parameter :: aa3 = 1.1_rkx
  !
  ! Number of gas taken into account by drydep scheme
  !
  integer(ik4), parameter :: ngasd = 31
  !
  ! threshold of rainfall intensity to activate water covered canopy option
  !
  real(rkx), parameter :: rainthr = 0.1_rkx
  !
  ! DATA section for the Zhang drydep scheme
  !
  ! The Zhang scheme uses its own LAI. First index is consistent with BATS
  ! LU type, 15 is for the number of month 12 + 3 for interp.
  ! The tables are supposed to be consistant with BATS types determined by
  ! ivegcov (be careffull with ethe ocean option)
  ! Rq: It would be better to have intercative directly : LAI and roughness
  ! from BATS and CLM
  ! Rq2: Since ivegcov is defined even when clm is activated, the dry dep
  ! scheme could in principle be used with CLM.
  ! BUT , there is also the option of activating CLM PFT level dydep scheme.
  !
  !NOTENOTENOTENOTENOTENOTENOTENOTEONOTENOTENOTENOTENOTENOTENOTENOTENOTENOTENOTE
  !
  ! NOTE : Now BATS, albeit not directly in the GLCC dataset but using the
  !        FUDGE, allows for 22 classes: urban and suburban (21 and 22) have
  !        been added. Here the scheme is inconsistent (20 classes only)
  !        Be aware of this.
  !
  !NOTENOTENOTENOTENOTENOTENOTENOTEONOTENOTENOTENOTENOTENOTENOTENOTENOTENOTENOTE
  !
  real(rkx) lai(22,15)
  real(rkx) z01(22) , z02(22)
  integer(ik4) :: kk

  data (lai(1,kk), kk = 1, 15)/                    &
           0.1_rkx , 0.1_rkx , 0.1_rkx , 0.5_rkx , 1.0_rkx , &
           2.0_rkx , 3.0_rkx , 3.5_rkx , 4.0_rkx , 0.1_rkx , &
           0.1_rkx , 0.1_rkx , 0.1_rkx , 0.1_rkx , 4.0_rkx /

  data (lai(2,kk), kk = 1, 15)/                    &
           1.0_rkx , 1.0_rkx , 1.0_rkx , 1.0_rkx , 1.0_rkx , &
           1.0_rkx , 1.0_rkx , 1.0_rkx , 1.0_rkx , 1.0_rkx , &
           1.0_rkx , 1.0_rkx , 1.0_rkx , 1.0_rkx , 1.0_rkx /

  data (lai(3,kk), kk = 1, 15)/                    &
           5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , &
           5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , &
           5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx /

  data (lai(4,kk), kk = 1, 15)/                    &
           0.1_rkx , 0.1_rkx , 0.5_rkx , 1.0_rkx , 2.0_rkx , &
           4.0_rkx , 5.0_rkx , 5.0_rkx , 4.0_rkx , 2.0_rkx , &
           1.0_rkx , 0.1_rkx , 0.1_rkx , 0.1_rkx , 5.0_rkx /

  data (lai(5,kk), kk = 1, 15)/                    &
           0.1_rkx , 0.1_rkx , 0.5_rkx , 1.0_rkx , 2.0_rkx , &
           4.0_rkx , 5.0_rkx , 5.0_rkx , 4.0_rkx , 2.0_rkx , &
           1.0_rkx , 0.1_rkx , 0.1_rkx , 0.1_rkx , 5.0_rkx /

  data (lai(6,kk), kk = 1, 15)/                    &
           6.0_rkx , 6.0_rkx , 6.0_rkx , 6.0_rkx , 6.0_rkx , &
           6.0_rkx , 6.0_rkx , 6.0_rkx , 6.0_rkx , 6.0_rkx , &
           6.0_rkx , 6.0_rkx , 6.0_rkx , 6.0_rkx , 6.0_rkx /

  data (lai(7,kk), kk = 1, 15)/                    &
           0.5_rkx , 0.5_rkx , 0.5_rkx , 0.5_rkx , 0.5_rkx , &
           0.5_rkx , 1.0_rkx , 2.0_rkx , 2.0_rkx , 1.5_rkx , &
           1.0_rkx , 1.0_rkx , 0.5_rkx , 0.5_rkx , 2.0_rkx  /

  data (lai(8,kk), kk = 1, 15)/                    &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx /

  data (lai(9,kk), kk = 1, 15)/                    &
           1.0_rkx , 1.0_rkx , 0.5_rkx , 0.1_rkx , 0.1_rkx , &
           0.1_rkx , 0.1_rkx , 1.0_rkx , 2.0_rkx , 1.5_rkx , &
           1.5_rkx , 1.0_rkx , 1.0_rkx , 0.1_rkx , 2.0_rkx /

  data (lai(10,kk), kk = 1, 15)/                   &
           1.0_rkx , 1.0_rkx , 1.0_rkx , 1.0_rkx , 1.0_rkx , &
           1.0_rkx , 1.0_rkx , 1.0_rkx , 1.0_rkx , 1.0_rkx , &
           1.0_rkx , 1.0_rkx , 1.0_rkx , 1.0_rkx , 1.0_rkx /

  data (lai(11,kk), kk = 1, 15)/                   &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx /

  data (lai(12,kk), kk = 1, 15)/                   &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx /

  data (lai(13,kk), kk = 1, 15)/                   &
           4.0_rkx , 4.0_rkx , 4.0_rkx , 4.0_rkx , 4.0_rkx , &
           4.0_rkx , 4.0_rkx , 4.0_rkx , 4.0_rkx , 4.0_rkx , &
           4.0_rkx , 4.0_rkx , 4.0_rkx , 4.0_rkx , 4.0_rkx /

  data (lai(14,kk), kk = 1, 15)/                   &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx /

  data (lai(15,kk), kk = 1, 15)/                   &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx /

  data (lai(16,kk), kk = 1, 15)/                   &
           3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx , &
           3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx , &
           3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx /

  data (lai(17,kk), kk = 1, 15)/                   &
           3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx , &
           3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx , &
           3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx /

  data (lai(18,kk), kk = 1, 15)/                  &
           3.0_rkx , 3.0_rkx , 3.0_rkx , 4.0_rkx , 4.5_rkx ,&
           5.0_rkx , 5.0_rkx , 5.0_rkx , 4.0_rkx , 3.0_rkx ,&
           3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx , 5.0_rkx /

  data (lai(19,kk), kk = 1, 15)/                  &
           3.0_rkx , 3.0_rkx , 3.0_rkx , 4.0_rkx , 4.5_rkx ,&
           5.0_rkx , 5.0_rkx , 5.0_rkx , 4.0_rkx , 3.0_rkx ,&
           3.0_rkx , 3.0_rkx , 3.0_rkx , 3.0_rkx , 5.0_rkx /

  data (lai(20,kk), kk = 1, 15)/                  &
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx ,&
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx ,&
           0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx /

  data  z01/0.02_rkx, 0.04_rkx, 0.90_rkx, 0.40_rkx, 0.40_rkx, 2.00_rkx,&
            0.02_rkx, 0.04_rkx, 0.03_rkx, 0.05_rkx, 0.04_rkx, 0.01_rkx,&
            0.10_rkx, 0.00_rkx, 0.00_rkx, 0.20_rkx, 0.20_rkx, 0.60_rkx,&
            0.60_rkx, 0.00_rkx, 0.00_rkx, 0.00_rkx/

  data  z02/0.10_rkx, 0.04_rkx, 0.90_rkx, 0.90_rkx, 1.00_rkx, 2.00_rkx,&
            0.10_rkx, 0.04_rkx, 0.03_rkx, 0.05_rkx, 0.04_rkx, 0.01_rkx,&
            0.10_rkx, 0.00_rkx, 0.00_rkx, 0.20_rkx, 0.20_rkx, 0.90_rkx,&
            0.90_rkx, 0.00_rkx, 0.00_rkx, 0.00_rkx/
  !
  ! Zhang stomatal resistance parameters
  !
  real(rkx) :: tmin(22) , tmax(22)
  real(rkx) :: rsminz(22) , brs(22)
  real(rkx) :: topt(22) , bvpd(22)
  real(rkx) :: psi1(22) , psi2(22)
  real(rkx) :: rac1(22) , rac2(22)
  real(rkx) :: rgo(22) , rcutdO(22)
  real(rkx) :: rcutwO(22) , rcutdS(22)
  real(rkx) :: rgs(22) , sdmax(22)
  real(rkx) :: mw(31) , rm(31)
  real(rkx) :: alphaz(31) , betaz(31)

  data tmin /  5.0_rkx,    5.0_rkx,   -5.0_rkx,   -5.0_rkx,    0.0_rkx, &
               0.0_rkx,    5.0_rkx, -999.0_rkx,   -5.0_rkx,    5.0_rkx, &
            -999.0_rkx, -999.0_rkx,    0.0_rkx, -999.0_rkx, -999.0_rkx, &
               0.0_rkx,    0.0_rkx,   -3.0_rkx,    0.0_rkx, -999.0_rkx, &
            -999.0_rkx,-999.0_rkx/

  data tmax / 45.0_rkx,   40.0_rkx,   40.0_rkx,   40.0_rkx,   45.0_rkx, &
              45.0_rkx,   45.0_rkx, -999.0_rkx,   40.0_rkx,   45.0_rkx, &
            -999.0_rkx, -999.0_rkx,   45.0_rkx, -999.0_rkx, -999.0_rkx, &
              45.0_rkx,   45.0_rkx,   42.0_rkx,   45.0_rkx, -999.0_rkx, &
            -999.0_rkx,-999.0_rkx/

  data rsminz / 120.0_rkx, 150.0_rkx, 250.0_rkx, 250.0_rkx, 150.0_rkx, &
                150.0_rkx, 100.0_rkx,-999.0_rkx, 150.0_rkx, 150.0_rkx, &
               -999.0_rkx,-999.0_rkx, 150.0_rkx,-999.0_rkx,-999.0_rkx, &
                150.0_rkx, 250.0_rkx, 150.0_rkx, 150.0_rkx,-999.0_rkx, &
               -999.0_rkx,-999.0_rkx/

  data brs / 40.0_rkx,  50.0_rkx,  44.0_rkx,  44.0_rkx,  43.0_rkx, &
             40.0_rkx,  20.0_rkx,-999.0_rkx,  25.0_rkx,  40.0_rkx, &
           -999.0_rkx,-999.0_rkx,  40.0_rkx,-999.0_rkx,-999.0_rkx, &
             40.0_rkx,  44.0_rkx,  44.0_rkx,  43.0_rkx,-999.0_rkx, &
           -999.0_rkx,-999.0_rkx/

  data topt / 27.0_rkx,  30.0_rkx,  15.0_rkx,  15.0_rkx,  27.0_rkx, &
              30.0_rkx,  25.0_rkx,-999.0_rkx,  20.0_rkx,  25.0_rkx, &
            -999.0_rkx,-999.0_rkx,  20.0_rkx,-999.0_rkx,-999.0_rkx, &
              30.0_rkx,  25.0_rkx,  21.0_rkx,  25.0_rkx,-999.0_rkx, &
            -999.0_rkx,-999.0_rkx/

  data bvpd /  0.00_rkx,   0.00_rkx,   0.31_rkx,   0.31_rkx,   0.36_rkx, &
               0.27_rkx,   0.00_rkx,-999.00_rkx,   0.24_rkx,   0.09_rkx, &
            -999.00_rkx,-999.00_rkx,   0.27_rkx,-999.00_rkx,-999.00_rkx, &
               0.27_rkx,   0.27_rkx,   0.34_rkx,   0.31_rkx,-999.00_rkx, &
            -999.0_rkx,-999.0_rkx/

  data psi1 / -1.5_rkx,  -1.5_rkx,  -2.0_rkx,  -2.0_rkx,  -1.9_rkx, &
              -1.0_rkx,  -1.5_rkx,-999.0_rkx,   0.0_rkx,  -1.5_rkx, &
            -999.0_rkx,-999.0_rkx,  -1.5_rkx,-999.0_rkx,-999.0_rkx, &
              -2.0_rkx,  -2.0_rkx,  -2.0_rkx,  -2.0_rkx,-999.0_rkx, &
            -999.0_rkx,-999.0_rkx/

  data psi2 / -2.5_rkx,  -2.5_rkx,  -2.5_rkx,  -2.5_rkx,  -2.5_rkx, &
              -5.0_rkx,  -2.5_rkx,-999.0_rkx,  -1.5_rkx,  -2.5_rkx, &
            -999.0_rkx,-999.0_rkx,  -2.5_rkx,-999.0_rkx,-999.0_rkx, &
              -4.0_rkx,  -3.5_rkx,  -2.5_rkx,  -3.0_rkx,-999.0_rkx, &
            -999.0_rkx,-999.0_rkx/

  data rac1 / 10.0_rkx, 20.0_rkx,100.0_rkx, 60.0_rkx, 40.0_rkx, &
             250.0_rkx, 10.0_rkx,  0.0_rkx, 40.0_rkx, 20.0_rkx, &
               0.0_rkx,  0.0_rkx, 20.0_rkx,  0.0_rkx,  0.0_rkx, &
              60.0_rkx, 40.0_rkx,100.0_rkx,100.0_rkx,  0.0_rkx, &
               0.0_rkx,  0.0_rkx/

  data rac2 / 40.0_rkx, 20.0_rkx,100.0_rkx,100.0_rkx, 40.0_rkx, &
             250.0_rkx, 40.0_rkx,  0.0_rkx, 40.0_rkx, 20.0_rkx, &
               0.0_rkx,  0.0_rkx, 20.0_rkx,  0.0_rkx,  0.0_rkx, &
              60.0_rkx, 40.0_rkx,100.0_rkx,100.0_rkx,  0.0_rkx, &
               0.0_rkx,  0.0_rkx/

  data rcutdO / 4000.0_rkx, 4000.0_rkx, 4000.0_rkx, 4000.0_rkx, 6000.0_rkx, &
                6000.0_rkx, 4000.0_rkx, -999.0_rkx, 8000.0_rkx, 4000.0_rkx, &
                -999.0_rkx, -999.0_rkx, 5000.0_rkx, -999.0_rkx, -999.0_rkx, &
                6000.0_rkx, 5000.0_rkx, 4000.0_rkx, 4000.0_rkx, -999.0_rkx, &
                -999.0_rkx, -999.0_rkx/

  data rcutwO / 200.0_rkx, 200.0_rkx, 200.0_rkx, 200.0_rkx, 400.0_rkx, &
                400.0_rkx, 200.0_rkx,-999.0_rkx, 400.0_rkx, 200.0_rkx, &
               -999.0_rkx,-999.0_rkx, 300.0_rkx,-999.0_rkx,-999.0_rkx, &
                400.0_rkx, 300.0_rkx, 200.0_rkx, 200.0_rkx,-999.0_rkx, &
               -999.0_rkx,-999.0_rkx/

  data rgO / 200.0_rkx,  200.0_rkx,  200.0_rkx,  200.0_rkx,  200.0_rkx, &
             200.0_rkx,  200.0_rkx,  500.0_rkx,  500.0_rkx,  500.0_rkx, &
             500.0_rkx, 2000.0_rkx,  500.0_rkx, 2000.0_rkx, 2000.0_rkx, &
             200.0_rkx,  200.0_rkx,  200.0_rkx,  200.0_rkx, 2000.0_rkx, &
               0.0_rkx,    0.0_rkx/

  data rcutds / 1500.0_rkx, 1000.0_rkx, 2000.0_rkx, 2000.0_rkx, 2500.0_rkx, &
                2500.0_rkx, 1000.0_rkx, -999.0_rkx, 2000.0_rkx, 2000.0_rkx, &
                -999.0_rkx, -999.0_rkx, 1500.0_rkx, -999.0_rkx, -999.0_rkx, &
                2000.0_rkx, 2000.0_rkx, 2500.0_rkx, 2500.0_rkx, -999.0_rkx, &
                   0.0_rkx,    0.0_rkx/

  data rgs / 200.0_rkx, 200.0_rkx, 200.0_rkx, 200.0_rkx, 200.0_rkx, &
             100.0_rkx, 200.0_rkx, 700.0_rkx, 300.0_rkx,  50.0_rkx, &
             700.0_rkx,  70.0_rkx,  50.0_rkx,  20.0_rkx,  20.0_rkx, &
             200.0_rkx, 200.0_rkx, 200.0_rkx, 200.0_rkx,  20.0_rkx, &
               0.0_rkx,   0.0_rkx/

  data sdmax /  10.0_rkx,  5.0_rkx, 200.0_rkx,   1.1_rkx, 200.0_rkx, &
               400.0_rkx, 20.0_rkx,   2.0_rkx,   2.0_rkx,  10.0_rkx, &
                 2.0_rkx,  1.0_rkx,  10.0_rkx,-999.0_rkx,-999.0_rkx, &
                50.0_rkx, 50.0_rkx, 200.0_rkx, 200.0_rkx,-999.0_rkx, &
                 0.0_rkx,  0.0_rkx/

  ! *****************************************************************
  ! * Gas Properties (Total 31 species)                          ****
  ! * Mesophyll resistance RM, scaling factors ALPHAZ and BETAZ, ****
  ! * molecular weight.                                          ****
  ! *****************************************************************
  ! parameters are given for the following species in the zhang scheme
  ! SO2    H2SO4   NO2    O3     H2O2   HONO  HNO3   HNO4  NH3
  ! PAN    PPN     APAN   MPAN
  ! HCHO   MCHO    PALD   C4A
  ! C7A    ACHO    MVK    MACR
  ! MGLY   MOH     ETOH   POH
  ! CRES   FORM    ACAC   ROOH
  ! ONIT   INIT

  data rm / 0.0_rkx ,   0.0_rkx ,   0.0_rkx ,   0.0_rkx ,    0.0_rkx , &
            0.0_rkx ,   0.0_rkx ,   0.0_rkx ,   0.0_rkx ,    0.0_rkx , &
            0.0_rkx ,   0.0_rkx ,   0.0_rkx ,   0.0_rkx ,  100.0_rkx , &
          100.0_rkx , 100.0_rkx , 100.0_rkx , 100.0_rkx ,    0.0_rkx , &
          100.0_rkx ,   0.0_rkx   , 0.0_rkx ,   0.0_rkx ,    0.0_rkx , &
            0.0_rkx ,   0.0_rkx   , 0.0_rkx ,   0.0_rkx ,  100.0_rkx , &
          100.0_rkx /

  data alphaz /  1.00_rkx , 1.00_rkx , 0.00_rkx , 0.00_rkx , 1.00_rkx , &
                10.00_rkx , 2.00_rkx , 5.00_rkx , 1.00_rkx , 0.00_rkx , &
                 0.00_rkx , 0.00_rkx , 0.00_rkx , 0.80_rkx , 0.00_rkx , &
                 0.00_rkx , 0.00_rkx , 0.00_rkx , 0.00_rkx , 0.00_rkx , &
                 0.00_rkx , 0.01_rkx , 0.60_rkx , 0.60_rkx , 0.40_rkx , &
                 0.01_rkx , 2.00_rkx , 1.50_rkx , 0.10_rkx , 0.00_rkx , &
                 0.00_rkx /

  data betaz /  0.00_rkx , 1.00_rkx , 0.80_rkx , 1.00_rkx , 1.00_rkx , &
               10.00_rkx , 2.00_rkx , 5.00_rkx , 0.00_rkx , 0.60_rkx , &
                0.60_rkx , 0.80_rkx , 0.30_rkx , 0.20_rkx , 0.05_rkx , &
                0.05_rkx , 0.05_rkx , 0.05_rkx , 0.05_rkx , 0.05_rkx , &
                0.05_rkx , 0.00_rkx , 0.10_rkx , 0.00_rkx , 0.00_rkx , &
                0.00_rkx , 0.00_rkx , 0.00_rkx , 0.80_rkx , 0.50_rkx , &
                0.50_rkx /


  data mw / 64.0_rkx ,  98.0_rkx ,  46.0_rkx ,  48.0_rkx ,  34.0_rkx , &
            63.0_rkx ,  47.0_rkx ,  79.0_rkx ,  17.0_rkx , 121.0_rkx , &
           135.0_rkx , 183.0_rkx , 147.0_rkx ,  30.0_rkx ,  44.0_rkx , &
            58.0_rkx ,  72.0_rkx , 128.0_rkx , 106.0_rkx ,  70.0_rkx , &
            70.0_rkx ,  72.0_rkx ,  32.0_rkx ,  46.0_rkx ,  60.0_rkx , &
           104.0_rkx ,  46.0_rkx ,  60.0_rkx ,  48.0_rkx ,  77.0_rkx , &
           147.0_rkx /

  contains

    subroutine drydep_aero(i,mbin,indsp,rhop,ivegcov,throw,roarow, &
                           ph,temp2,sutemp,srad,rh10,      &
                           wind10,zeff,beffdiam,pdepv,ddepv,ustar,ra)
      implicit none
      integer(ik4) , intent(in) :: i , mbin
      integer(ik4) , intent(in) , dimension(mbin) :: indsp
      integer(ik4) , intent(in) , dimension(jci1:jci2) :: ivegcov
      real(rkx) , dimension(jci1:jci2) , intent(in) :: rh10 , &
                       srad , sutemp , temp2 , wind10 , zeff
      real(rkx) , dimension(jci1:jci2,kz) , intent(in) :: ph , roarow , throw
      real(rkx) , dimension(jci1:jci2) , intent(in) :: ra , ustar
      real(rkx) , dimension(mbin) , intent(in) :: beffdiam
      real(rkx) , intent(in) :: rhop

      ! output table to be passed out. Care dimension is ntr

      real(rkx) , intent(out) , dimension(jci1:jci2,kz,ntr) :: pdepv
      real(rkx) , intent(out) , dimension(jci1:jci2,ntr) :: ddepv

      real(rkx) :: amfp , amob , eb , eim , ein , frx1 , w1 , w2
      real(rkx) :: pre , prii , priiv , r1 , st , rhsize , pdiff
      real(rkx) , dimension(jci1:jci2,kz) :: amu
      real(rkx) , dimension(jci1:jci2) :: anu , schm
      real(rkx) , dimension(jci1:jci2,kz,mbin) :: cfac , pdepvsub , taurel
      real(rkx) , dimension(jci1:jci2,mbin) :: rs
      real(rkx), dimension(jci1:jci2,2:kz) :: wk, settend
      real(rkx) , dimension(mbin) :: avesize
      integer(ik4) :: j , k , lcov , l , n , ib
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'drydep_aero'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      ! here avesize is a RADIUS of the particle bin in m :
      ! calculated here from bin effective diameter in micrometer
      do n = 1 , mbin
        avesize(n) = (beffdiam(n) * d_half) * 1.e-6_rkx
      end do
      ! ********************************************************
      ! *   aerosize - dry radius                    !!     ****
      ! *   rhop  - density for each aerosol type           ****
      ! ********************************************************
      do n = 1 , mbin
        do k = 1 , kz
          do j = jci1 , jci2
            !
            ! ********************************************************
            ! *  aerosol gravitational settling velocity          ****
            ! *  and diffusion coefficient                        ****
            ! *                                                   ****
            ! * air's dynamic viscosity                           ****
            ! * Sutherland Equation
            ! ********************************************************
            !
            amu(j,k) = (a1*(throw(j,k)**a2))/(throw(j,k)+a3)
            ! mid layer pressure in [pascal].
            pre = ph(j,k)
            !
            ! ********************************************************
            ! * mean molecular free path.                         ****
            ! *     K.V. Beard [1976], J Atm. Sci., 33            ****
            ! ********************************************************
            !
            amfp = c1*(amu(j,k)/c2)*(c3/pre)*sqrt(throw(j,k)/c4)
            prii = 2.0_rkx/9.0_rkx*egrav/amu(j,k)
            priiv = prii*(rhop-roarow(j,k))
            !
            ! ********************************************************
            ! * Cunningham slip correction factor and             ****
            ! * relaxation time = vg/grav.                        ****
            ! ********************************************************
            !
            cfac(j,k,n) = d_one + amfp/avesize(n) * &
                         (aa1+aa2*exp(-aa3*avesize(n)/amfp))
            taurel(j,k,n) = priiv*avesize(n)**2*cfac(j,k,n) * regrav
            !
            ! ********************************************************
            ! * stokes friction                                  *****
            ! ! pdepvsub(j,k,n) ' sellting dep. velocity = '
            ! ********************************************************
            !
            pdepvsub(j,k,n) = taurel(j,k,n)*egrav
          end do
        end do
      end do
      !
      ! *****************************************************
      ! * the schmidt number is the ratio of the         ****
      ! * kinematic viscosity of air to the particle     ****
      ! * brownian diffusivity ===> sc=v/d               ****
      ! *****************************************************
      !
      do n = 1 , mbin
        do j = jci1 , jci2
          !
          ! *****************************************************
          ! * for now we will not consider the humidity      ****
          ! * impact so we will set the variable frx1=1.0    ****
          ! * i.e. only dry particles                        ****
          ! *****************************************************
          !
          frx1 = 1.0_rkx
          rhsize = avesize(n)*frx1 ! still a radius
          anu(j) = amu(j,kz)/roarow(j,kz)
          amob = 6.0_rkx*mathpi*amu(j,kz)*rhsize/cfac(j,kz,n)
          pdiff = boltzk*throw(j,kz)/amob
          schm(j) = anu(j)/pdiff
          !
          ! ******************************************************
          ! * for brownian diffusion, there is evidence that  ****
          ! * its fromula depend on schmidt number as :       ****
          ! * eb= schm x c^gama                               ****
          ! * where gama is efficiency factor and its value   ****
          ! * between 1/2 and 2/3 with larger values for      ****
          ! * rougher surfaces                                ****
          ! * ****************************************************
          !
        end do
        do j = jci1 , jci2
            !
            ! find the right table index for the cell cover ( ocean
            ! and lake are 0 in the ivegcov and 14-15 in the table )
            !
            if ( ivegcov(j) == 0 ) then
              lcov = 14
            else if ( ivegcov(j) > 20 ) then
              lcov = 20
            else
              lcov = ivegcov(j)
            end if
            !
            ! ******************************************************
            ! * the parameter governing impaction processes is *****
            ! * the stokes number,st, which has the form of    *****
            ! * 1) st = vg x u* /g a for vegetated surefaces   *****
            ! *    (slinn, 1982)                               *****
            ! * 2) st = vg x u*2/anu for smothed surfaces or   *****
            ! *    surfaces with bluff roughness elements      *****
            ! *    (giorgi,1988)                               *****
            ! ******************************************************
            !
            ! Graziano - 2018-02-09 - updated Stokes number computation
            ! Before, only the formulation for smooth surfaces was used
            if ( ast(lcov) > d_zero ) then
              st = taurel(j,kz,n)*ustar(j)*regrav/ast(lcov)
              eb = schm(j)**(-agam(lcov))
            else
              st = taurel(j,kz,n)*ustar(j)*ustar(j)/anu(j)
              !eb = schm(j)**(-twot)
              eb = schm(j)**(-d_half)
            end if
            eim = (st/(st+aest(lcov)))**2
            eim = max(min(eim,0.6_rkx),1.0e-08_rkx)
            if ( arye(lcov) > 0.001_rkx ) then
              ein = d_two*((1000.0_rkx*avesize(n))/arye(lcov))**1.5_rkx
              !ein = d_two*((1000.0_rkx*avesize(n))/arye(lcov))**2
            else
              ein = 1.0e-08_rkx
            end if
            ein = max(min(ein,0.5_rkx),1.0e-08_rkx)
            !
            ! *****************************************************
            ! * partickes larger than 5 micro may rebounded   *****
            ! * after hitting a surface,this process may be   *****
            ! * included by modifying the total collection    *****
            ! * by the factor of r1, which represents the     *****
            ! * fraction of particles sticking to the surface *****
            ! * slinn (1982) suggested the following:         *****
            ! * r = exp (- st^0.2)                            *****
            ! *****************************************************
            r1 = exp(-sqrt(st))
            if ( r1 < 0.4_rkx ) r1 = 0.4_rkx
            ! ***************************************************
            ! * calculation of rs: the surface resistance   *****
            ! * which depends on the collection efficiency  *****
            ! * of the surface and is determined by the     *****
            ! * various deposition processes                *****
            ! ***************************************************
            rs(j,n) = 3.0_rkx*ustar(j)*(eb+eim+ein)*r1
            rs(j,n) = max(1.0e-5_rkx,min(1.e5_rkx,rs(j,n)))
            rs(j,n) = d_one/rs(j,n)
        end do
      end do

      ! average settling and deposition velocities on bin
      ! care we use pdepv and ddpv table that are dimensionned to ntr
      ! and not mbin !
      do ib = 1 , mbin
        ! there isw no sub-bin anymore / we consider directly effective radius
        pdepv(:,:,indsp(ib)) = 0.0_rkx
        ddepv(:,indsp(ib))   = 0.0_rkx
        do k = 1 , kz
          do j = jci1 , jci2
            pdepv(j,k,indsp(ib)) = pdepvsub(j,k,ib)
          end do
        end do
        do j = jci1 , jci2
            ! agregate the dry deposition velocity, remember one cover per grid
            ! cell for now
            ! the dry deposition velocity must account also for the
            ! settling velocity at kz
            ! simple form now add the vs
            ddepv(j,indsp(ib)) = 1.0_rkx/(ra(j)+rs(j,ib)) + &
                       pdepvsub(j,kz,ib)
        end do
      end do
      !
      ! Finally update the emission and settling tendencies for
      ! dust and sea salt
      !
      if ( idynamic == 3 ) then
        do ib = 1 , mbin
          ! deposition, remember chiten must be consistent with chemt
          do k = 2 , kz
            do j = jci1 , jci2
              if ( chemt(j,i,k-1,indsp(ib)) > mintr ) then
                w1 = cfmz(j,i,k)/cfmz(j,i,k-1)
                w2 = d_two - w1
                wk(j,k) = d_half * (w2 * chemt(j,i,k,indsp(ib)) + &
                                    w1 * chemt(j,i,k-1,indsp(ib))) * rdt
              else
                wk(j,k) = d_zero
              end if
            end do
          end do
          do j = jci1 , jci2
            do k = 2 , kz - 1
              ! do not apply to the first level
              !settend(j,k) = (wk(j,k+1)*pdepv(j,k+1,indsp(ib)) - &
              !                wk(j,k)*pdepv(j,k,indsp(ib))) / cdzq(j,i,k)
              ! use exponential form for stability
              settend(j,k) =  wk(j,k+1) * &
                  (d_one - exp(-pdepvsub(j,k+1,ib)/cdzq(j,i,k)*dt)) - &
                              wk(j,k)   * &
                  (d_one - exp(-pdepvsub(j,k,ib)/cdzq(j,i,k)*dt))
              chiten(j,i,k,indsp(ib)) = chiten(j,i,k,indsp(ib)) - settend(j,k)
              if ( ichdiag > 0 ) then
                cseddpdiag(j,i,k,indsp(ib)) = cseddpdiag(j,i,k,indsp(ib)) - &
                                              settend(j,k) * cfdout
              end if
            end do
            !
            ! option 1 : calculate the tend as flux divergence
            ! at first level include surface drydep velocity to calculate the
            ! divergence
            !
            ! option 2 : the dry deposition is accounted for in the BL scheme
            ! we just pass the surface flux to the pbl interface (actually the
            ! net surface flux, cf also emission module)
            !
            if ( ichdrdepo == 1 ) then
              !settend(j,kz) = (chib(j,i,kz,indsp(ib))  * ddepv(j,indsp(ib))-  &
              !                 wk(j,kz)*pdepv(j,kz,indsp(ib))) / cdzq(j,i,kz)
              ! use exponential form for stability
              settend(j,kz) = max(chemt(j,i,kz,indsp(ib)),d_zero)*rdt * &
                  (d_one - exp(-ddepv(j,indsp(ib))/cdzq(j,i,kz)*dt)) - &
                                wk(j,kz) * &
                  (d_one - exp(-pdepvsub(j,kz,ib)/cdzq(j,i,kz)*dt))
              chiten(j,i,kz,indsp(ib)) = chiten(j,i,kz,indsp(ib)) - &
                  settend(j,kz)
              ! save the dry deposition flux for coupling with
              ! landsurface scheme (Kg.m2.s-1)
              ! consider ddflux = Cav . Vd where Cav would be the average
              ! concentration within the time step Cav = 0.5 (C + (C+deltaC))
              ! care  chib and settend have to be corrected for pressure
              ! cdrydepflux is a time accumulated array set to zero when surface
              ! scheme is called (cf atm to surf interface)
              cdrydepflx(j,i,indsp(ib)) = cdrydepflx(j,i,indsp(ib)) + &
                (chemt(j,i,kz,indsp(ib))  - settend(j,kz)*dt/d_two) * &
                   crhob3d(j,i,kz) * ddepv(j,indsp(ib))

              !diagnostic for settling and drydeposition removal
              if ( ichdiag > 0 ) then
                cseddpdiag(j,i,kz,indsp(ib)) = cseddpdiag(j,i,kz,indsp(ib)) - &
                                               settend(j,kz) * cfdout
              end if

              ! accumulated diagnostic for dry deposition flux
              ! average (in kg .m2.s-1)
              ! remdrd(j,i,indsp(ib)) = remdrd(j,i,indsp(ib)) + &
              !                         cdrydepflx(j,i,indsp(ib)) * cfdout
              ! bugfix: remdrd is accumulated and averaged at a different
              ! frequency than drydepflx accum !
              ! do not use drydepflx in the formula
              remdrd(j,i,indsp(ib)) = remdrd(j,i,indsp(ib)) + &
                       (chemt(j,i,kz,indsp(ib)) - settend(j,kz)*dt/d_two) * &
                        crhob3d(j,i,kz)*ddepv(j,indsp(ib)) * cfdout
              ! alternative formulation using tendency/flux relationship
              ! remdrd(j,i,indsp(ib)) = remdrd(j,i,indsp(ib)) + &
              !            chemt(j,i,kz,indsp(ib)) * &
              !         (d_one - exp(-ddepv(j,indsp(ib)) / &
              !                 cdzq(j,i,kz)*dt ))*rdt  &
              !           * crhob3d(j,i,kz) / cdzq(j,i,kz) * cfdout

              ! no net flux is passed to BL schemes in this case
              chifxuw(j,i,indsp(ib)) = d_zero
              drydepv(j,i,indsp(ib)) = d_zero
            else if ( ichdrdepo == 2 ) then
              !
              ! add the dry deposition term to the net emision/deposition flux
              ! for the BL scheme !
              ! flux
              chifxuw(j,i,indsp(ib)) = chifxuw(j,i,indsp(ib)) - &
                  max(chemt(j,i,kz,indsp(ib)),d_zero)*rdt * &
                  (d_one - exp(-ddepv(j,indsp(ib))/cdzq(j,i,kz)*dt)) + &
                       wk(j,kz) * &
                   (d_one - exp(-pdepvsub(j,kz,ib)/cdzq(j,i,kz)*dt))
              drydepv(j,i,indsp(ib)) = ddepv(j,indsp(ib))
            end if
            !
            ! dry dep velocity diagnostic in m.s-1  ( + drydep v. include
            ! also settling , accumulated between two outputs time step)
            ddv_out(j,i,indsp(ib)) = ddv_out(j,i,indsp(ib)) + &
                   ddepv(j,indsp(ib))
          end do
        end do
      else
        do ib = 1 , mbin
          ! deposition, remember chiten must be normalised by psb and
          ! consistent with chib
          do k = 2 , kz
            do j = jci1 , jci2
              if ( chib(j,i,k-1,indsp(ib)) > mintr * cpsb(j,i) ) then
                wk(j,k) = (twt(k,1)*chib(j,i,k,indsp(ib)) + &
                           twt(k,2)*chib(j,i,k-1,indsp(ib)))*rdt
              else
                wk(j,k) = d_zero
              end if
            end do
          end do
          do j = jci1 , jci2
            do k = 2 , kz - 1
              ! do not apply to the first level
              !settend(j,k) = (wk(j,k+1)*pdepv(j,k+1,indsp(ib)) - &
              !                wk(j,k)*pdepv(j,k,indsp(ib))) / cdzq(j,i,k)
              ! use exponential form for stability
              settend(j,k) =  wk(j,k+1) * &
                  (d_one - exp(-pdepvsub(j,k+1,ib)/cdzq(j,i,k)*dt)) - &
                              wk(j,k)   * &
                  (d_one - exp(-pdepvsub(j,k,ib)/cdzq(j,i,k)*dt))
              chiten(j,i,k,indsp(ib)) = chiten(j,i,k,indsp(ib)) - settend(j,k)
              if ( ichdiag > 0 ) then
                cseddpdiag(j,i,k,indsp(ib)) = cseddpdiag(j,i,k,indsp(ib)) - &
                                              settend(j,k) * cfdout
              end if
            end do
            !
            ! option 1 : calculate the tend as flux divergence
            ! at first level include surface drydep velocity to calculate the
            ! divergence
            !
            ! option 2 : the dry deposition is accounted for in the BL scheme
            ! we just pass the surface flux to the pbl interface (actually the
            ! net surface flux, cf also emission module)
            !
            if ( ichdrdepo == 1 ) then
              !settend(j,kz) = (chib(j,i,kz,indsp(ib))  * ddepv(j,indsp(ib))-  &
              !                 wk(j,kz)*pdepv(j,kz,indsp(ib))) / cdzq(j,i,kz)
              ! use exponential form for stability
              settend(j,kz) = max(chib(j,i,kz,indsp(ib)),d_zero)*rdt * &
                  (d_one - exp(-ddepv(j,indsp(ib))/cdzq(j,i,kz)*dt)) - &
                                wk(j,kz) * &
                  (d_one - exp(-pdepvsub(j,kz,ib)/cdzq(j,i,kz)*dt))
              chiten(j,i,kz,indsp(ib)) = chiten(j,i,kz,indsp(ib)) - &
                  settend(j,kz)
              ! save the dry deposition flux for coupling with
              ! landsurface scheme (Kg.m2.s-1)
              ! consider ddflux = Cav . Vd where Cav would be the average
              ! concentration within the time step Cav = 0.5 (C + (C+deltaC))
              ! care  chib and settend have to be corrected for pressure
              ! cdrydepflux is a time accumulated array set to zero when surface
              ! scheme is called (cf atm to surf interface)
              cdrydepflx(j,i,indsp(ib)) = cdrydepflx(j,i,indsp(ib)) + &
                (chib(j,i,kz,indsp(ib))  - settend(j,kz)*dt/d_two) / &
                   cpsb(j,i) * crhob3d(j,i,kz) * ddepv(j,indsp(ib))

              !diagnostic for settling and drydeposition removal
              if ( ichdiag > 0 ) then
                cseddpdiag(j,i,kz,indsp(ib)) = cseddpdiag(j,i,kz,indsp(ib)) - &
                                               settend(j,kz) * cfdout
              end if

              ! accumulated diagnostic for dry deposition flux
              ! average (in kg .m2.s-1)
              ! remdrd(j,i,indsp(ib)) = remdrd(j,i,indsp(ib)) + &
              !                         cdrydepflx(j,i,indsp(ib)) * cfdout
              ! bugfix: remdrd is accumulated and averaged at a different
              ! frequency than drydepflx accum !
              ! do not use drydepflx in the formula
              remdrd(j,i,indsp(ib)) = remdrd(j,i,indsp(ib)) + &
                       (chib(j,i,kz,indsp(ib)) - settend(j,kz)*dt/d_two) / &
                        cpsb(j,i) * crhob3d(j,i,kz)*ddepv(j,indsp(ib)) * cfdout
              ! alternative formulation using tendency/flux relationship
              ! remdrd(j,i,indsp(ib)) = remdrd(j,i,indsp(ib)) + &
              !            chib3d(j,i,kz,indsp(ib)) * &
              !         (d_one - exp(-ddepv(j,indsp(ib)) / &
              !                 cdzq(j,i,kz)*dt ))*rdt  &
              !           * crhob3d(j,i,kz) / cdzq(j,i,kz) * cfdout

              ! no net flux is passed to BL schemes in this case
              chifxuw(j,i,indsp(ib)) = d_zero
              drydepv(j,i,indsp(ib)) = d_zero
            else if ( ichdrdepo == 2 ) then
              !
              ! add the dry deposition term to the net emision/deposition flux
              ! for the BL scheme !
              ! flux
              chifxuw(j,i,indsp(ib)) = chifxuw(j,i,indsp(ib)) - &
                  max(chib(j,i,kz,indsp(ib)),d_zero)*rdt * &
                  (d_one - exp(-ddepv(j,indsp(ib))/cdzq(j,i,kz)*dt)) + &
                       wk(j,kz) * &
                   (d_one - exp(-pdepvsub(j,kz,ib)/cdzq(j,i,kz)*dt))
              chifxuw(j,i,indsp(ib)) = chifxuw(j,i,indsp(ib))/cpsb(j,i)
              drydepv(j,i,indsp(ib)) = ddepv(j,indsp(ib))
            end if
            !
            ! dry dep velocity diagnostic in m.s-1  ( + drydep v. include
            ! also settling , accumulated between two outputs time step)
            ddv_out(j,i,indsp(ib)) = ddv_out(j,i,indsp(ib)) + &
                   ddepv(j,indsp(ib))
          end do
        end do
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine drydep_aero

    subroutine drydep_gas(i,lmonth,lday,ivegcov,rh10,srad,tsurf, &
                          prec,temp10,xlai,ustar,resa)
      use mod_che_indices
      implicit none
      integer(ik4) , intent(in) :: i
      integer(ik4), intent(in) :: lmonth , lday
      integer(ik4) , intent(in) , dimension(jci1:jci2) :: ivegcov
      real(rkx) , intent(in) , dimension(jci1:jci2) :: rh10 , srad , tsurf , &
                                            prec, temp10, xlai 
      real(rkx) , dimension(jci1:jci2) , intent(in) :: ustar , resa
      real(rkx),  dimension(jci1:jci2,ntr) :: drydepvg

      integer(ik4) :: n , j , im , l , lcov
      real(rkx) , dimension(ngasd,jci1:jci2) :: resb, resc
      real(rkx) , dimension(ngasd,jci1:jci2) :: vdg
      real(rkx) , dimension(jci1:jci2) :: icz , ddrem
      real(rkx) , dimension(jci1:jci2) :: lai_f , laimin , laimax , snow
      real(rkx) :: kd , kav , rdz
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'drydep_gas'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      ! Different options for LAI and roughness
      ! for the moment read from

      do j = jci1 , jci2
        if ( ivegcov(j) == 0 ) then
          lcov = 14
        else if ( ivegcov(j) > 20 ) then
          lcov = 20
        else
          lcov = ivegcov(j)
        end if

!        im = lmonth - 1
!        if ( lmonth == 1 ) im = 12
!        if (lday <= 15 ) then
!          lai_f(j) = lai(lcov,im) + (lai(lcov,lmonth) - &
!                     lai(lcov,im))/30._rkx * real(15 + lday,rkx)
!        else
!          lai_f(j) = lai(lcov,lmonth) + (lai(lcov,lmonth+1) - &
!                     lai(lcov,lmonth))/30._rkx * real(lday - 15,rkx)
!        end if
!        if ( lai_f(j) < d_zero) lai_f(j) = d_zero
        laimin(j) = lai(lcov,14)
        laimax(j) = lai(lcov,15)
! FAB use interactive LAI by default,laimin and laimax are kept for further scaling. 
        lai_f(j)  = xlai(j)  
      end do
     

      snow(:) = d_zero
      icz(:) = czen(:,i)
      call stomtresis(lai_f,laimin,laimax,ivegcov,ngasd,ustar,prec,snow,srad, &
                      tsurf,temp10,rh10,icz,resc,resb)
      ! now calculate the dry deposition velocities and select it
      ! according to the gasphase mechanism
      ! vdg in m.s-1
      vdg(:,:) = d_zero
      do j = jci1 , jci2
          do n = 1 , ngasd
            vdg(n,j) = d_one/(resa(j)+resb(n,j)+resc(n,j))
          end do
      end do
      ! this part depends on the chem mechanism
      ! for CBMZ , we can certainly improve this.
      drydepvg = d_zero
      drydepvg(jci1:jci2,iso2)  =  vdg(1,jci1:jci2)
      ! SO2 deposition is used in SULF , AERO and CBMZ simulations
      if ( igaschem > 0 ) then
        drydepvg(jci1:jci2,ino2)  =  vdg(3,jci1:jci2)!*0.5
        drydepvg(jci1:jci2,io3)   =  vdg(4,jci1:jci2)!*0.5
        drydepvg(jci1:jci2,ih2o2) =  vdg(5,jci1:jci2)!*0.5
        drydepvg(jci1:jci2,ihno3) =  vdg(6,jci1:jci2)!*0.5
        drydepvg(jci1:jci2,inh3)  =  vdg(9,jci1:jci2)!*0.5
        drydepvg(jci1:jci2,ipan)  =  vdg(10,jci1:jci2)!*0.5
        drydepvg(jci1:jci2,ihcho) =  vdg(14,jci1:jci2)!*0.5
        drydepvg(jci1:jci2,iald2) =  vdg(15,jci1:jci2)!*0.5
        drydepvg(jci1:jci2,ich3oh)  =  vdg(23,jci1:jci2)!*0.5

        !FAB try out wesely from CLM 
      !  drydepvg(jci1:jci2,io3) = cddepv_clm(jci1:jci2,i,4)  
      end if

      ! Finally : gas phase dry dep tendency calculation
      if ( ichdrdepo == 1 ) then
        do j = jci1 , jci2
          rdz = d_one / cdzq(j,i,kz)
          do n = 1 , ntr
            kd = drydepvg(j,n) * rdz !Kd removal rate in s-1
            if ( idynamic == 3 ) then
              kav = max(chemt(j,i,kz,n),d_zero)*rdt
            else
              kav = max(chib3d(j,i,kz,n),d_zero)*rdt
            end if
            if ( kd*dt < 25.0_rkx ) then
              ! dry dep removal tendency (+)
              ddrem(j) = kav * (d_one-exp(-kd*dt))
            else
              ddrem(j) = d_zero
            end if
            ! update chiten
            chiten(j,i,kz,n) = chiten(j,i,kz,n) - ddrem(j)
            ! diag dry dep tendency
            if ( ichdiag > 0 ) then
               cseddpdiag(j,i,kz,n) = cseddpdiag(j,i,kz,n) - &
                                               ddrem(j) *  cfdout
            end if
            ! drydep flux diagnostic (accumulated between two outputs time
            ! step) ! flux is in kg/m2/s-1.
            remdrd(j,i,n) = remdrd(j,i,n) + max(chemt(j,i,kz,n),d_zero) &
                                          * crhob3d(j,i,kz)*drydepvg(j,n)* cfdout
            ! dry dep velocity diagnostic in m.s-1
            ! (accumulated between two outputs time step)
            drydepv(j,i,n) = d_zero
            ddv_out(j,i,n) = ddv_out(j,i,n) + drydepvg(j,n)
          end do
        end do
      else if ( ichdrdepo == 2 ) then
        do j = jci1 , jci2
          do n = 1 , ntr
            if ( idynamic == 3 ) then
              chifxuw(j,i,n) = chifxuw(j,i,n) - chemt(j,i,kz,n) * drydepvg(j,n)
            else
              chifxuw(j,i,n) = chifxuw(j,i,n) - chib3d(j,i,kz,n)* drydepvg(j,n)
            end if
            ! dry dep velocity diagnostic in m.s-1
            ! (accumulated between two outputs time step)
            drydepv(j,i,n) =  drydepvg(j,n)
            ddv_out(j,i,n) =  ddv_out(j,i,n) + drydepvg(j,n)
          end do
        end do
      end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine drydep_gas

    subroutine aerodyresis(zeff,wind10,temp2,sutemp,rh10,srad,ivegcov,ustar,ra)
      implicit none
      integer(ik4) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: ivegcov
      real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: temp2
      real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: wind10
      real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: rh10
      real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: sutemp
      real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: srad
      real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: zeff
      real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(out) :: ustar
      real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(out) :: ra
      integer(ik4) :: i , j , l
      real(rkx) :: vp , tsv
      real(rkx) :: z , zl , ww
      real(rkx) :: ptemp2 , es , qs
      real(rkx) :: wvpm , vptemp , tsw , mol
      real(rkx) :: z0water , dthv , cun , zdl
      real(rkx) :: psiu , psit , x , y
      real(rkx) :: thstar , rib , dtemp , tbar
      real(rkx) :: ustarsq , utstar , kui
      real(rkx) :: ratioz , logratio , asq
      real(rkx) :: aa , cm , ch , fm , fh , zz0
      real(rkx) , parameter :: z10 = 10.0_rkx
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'aerodyresis'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif

      ! ****************************************************
      ! * ra : is the aerodynamic resistance above the  ****
      ! *      canopy and it is function in u* and      ****
      ! *      z0: roughness length and the stability   ****
      ! *      function                                 ****
      ! * mol  - monin obukhov length (m) - calculated  ****
      ! *           for each land use category          ****
      ! * ptemp2 -potential temperature at z2  (deg. k) ****
      ! * temp2 - temperature at 10m. (deg k)           ****
      ! * z10   - 10 m.                                 ****
      ! * sutemp -surface temperature (deg k)           ****
      ! * srad   -solar irradiance at the ground(w/m**2)****
      ! * rh10  - relative humidity of air at 10m.      ****
      ! *           (0.0-1.0)                           ****
      ! * stdpmb - sea level pressure (mb)              ****
      ! ****************************************************
      do i = ici1 , ici2
        do j = jci1 , jci2
            ww = max(wind10(j,i),1.0_rkx)
            zz0 = zeff(j,i)
            ! ***************************************************************
            ! * potential temperature at z2  (deg. k)
            ! ***************************************************************
            ptemp2 = temp2(j,i) + z10*0.0098_rkx
            ! ***************************************************************
            ! * for calculations over water compute values of critical
            ! * profile variables: l and ustar
            ! *           ******begin for water***
            ! ***************************************************************
            if ( ivegcov(j,i) == 0 ) then
              ! **************************************************************
              ! * vp  - vapour pressure at z2
              ! * wvpm- water vapour mixing ratio at  z2
              ! * vptemp- virtual potential temperature at z2
              ! **************************************************************
              es = 6.108_rkx * &
                exp(17.27_rkx*(temp2(j,i)-tzero)/(temp2(j,i)-35.86_rkx))
              vp = rh10(j,i)*es
              wvpm = ep2*vp/(stdpmb-vp)
              vptemp = ptemp2*(1.0_rkx+0.61_rkx*wvpm)
              ! **************************************************************
              ! * assume rh10 at water surface is 100%
              ! *   vp = es(tsw-tzero) !sat. vap press at surface
              ! *   saturated vapour pressure at surface
              ! *   saturated mixing ratio at surface
              ! *   tsv - virtual potential temperature at surface (deg. k)
              ! **************************************************************
              tsw = sutemp(j,i)
              vp = 6.108_rkx*exp(17.27_rkx*(tsw-tzero)/(tsw-35.86_rkx))
              qs = ep2*vp/(stdpmb-vp)
              tsv = tsw*(1.0_rkx+0.61_rkx*qs)
              z0water = 1.0e-4_rkx
              ! **************************************************************
              ! * scalet  :  not required if  z2 = 10m
              ! **************************************************************
              dthv = (vptemp-tsv)
              ! **************************************************************
              ! * calculate drag coefficient cun with neutral condition
              ! * assumption  garratt (1977)
              ! **************************************************************
              cun = 7.5e-4_rkx + 6.7e-5_rkx*ww
              mol = 9999.0_rkx
              if ( abs(dthv) > 1.0e-6_rkx ) then
                mol = vptemp*cun**1.5_rkx*ww**2/(5.096e-3_rkx*dthv)
              end if
              if ( mol > 0.0_rkx  .and. mol < 5.0_rkx ) mol =  5.0_rkx
              if ( mol > -5.0_rkx .and. mol < 0.0_rkx ) mol = -5.0_rkx
              zdl = z10/mol
              if ( zdl < 0.0_rkx ) then
                ! **************************************************************
                ! * wind speed
                ! **************************************************************
                x = (1.0_rkx-15.0_rkx*zdl)**0.25_rkx
                psiu = 2.0_rkx*log(0.5_rkx*(1.0_rkx+x)) + &
                               log(0.5_rkx*(1.0_rkx+x*x)) - &
                       2.0_rkx*atan(x) + 0.5_rkx*mathpi
                ! **************************************************************
                ! * pot temp
                ! **************************************************************
                y = sqrt(1.0_rkx-9.0_rkx*zdl)
                psit = 2.0_rkx*0.74_rkx*log((1.0_rkx+y)/2.0_rkx)
              else
                psiu = -4.7_rkx*zdl
                psit = psiu
              end if
              z0water = 0.000002_rkx*ww**2.5_rkx
              ustar(j,i) = vonkar*ww/(log(z10/z0water)-psiu)
              thstar = vonkar*(ptemp2-sutemp(j,i)) / &
                       (0.74_rkx*log(z10/z0water)-psit)
              zz0 = z0water
            else
              ! **************************************************************
              ! * compute ustar and l for land use categories other than
              ! * water use louis method. !pkk 7/16/85, find bulk
              ! * richardson number.
              ! **************************************************************
              rib = egrav*z10*(ptemp2-sutemp(j,i))/(sutemp(j,i)*ww**2)
              ! ***************************************************************
              ! * ensure that conditions over land are never stable when
              ! * there is incoming solar radiation
              ! ***************************************************************
              if ( srad(j,i) > 0.0_rkx .and. rib > 0.0_rkx ) rib = 1.e-15_rkx
              dtemp = ptemp2 - sutemp(j,i)
              if ( abs(dtemp) < 1.e-10_rkx ) dtemp = sign(1.e-10_rkx,dtemp)
              tbar = 0.5_rkx*(ptemp2+sutemp(j,i))
              ratioz = z10/zz0
              logratio = log(ratioz)
              asq = 0.16_rkx/(logratio**2)
              if ( rib <= 0.0_rkx ) then
                aa = asq*9.4_rkx*sqrt(ratioz)
                cm = 7.4_rkx*aa
                ch = 5.3_rkx*aa
                fm = 1.0_rkx - (9.4_rkx*rib/(1.0_rkx+cm*sqrt(abs(rib))))
                fh = 1.0_rkx - (9.4_rkx*rib/(1.0_rkx+ch*sqrt(abs(rib))))
              else
                fm = 1.0_rkx/((1.0_rkx+4.7_rkx*rib)**2)
                fh = fm
              end if
              ustarsq = asq*ww**2*fm
              utstar = asq*ww*dtemp*fh/0.74_rkx
              ustar(j,i) = sqrt(ustarsq)
              thstar = utstar/ustar(j,i)
              mol = tbar*ustarsq/(vonkar*egrav*thstar)
            end if

            kui = 1.0_rkx/(vonkar*ustar(j,i))

            ! **************************************************************
            ! * compute the values of  ra                            *******
            ! **************************************************************
            z = z10
            zl = z/mol
            if ( zl >= 0.0_rkx ) then
              ra(j,i) = kui*(0.74_rkx*log(z/zz0)+4.7_rkx*zl)
            else
              ra(j,i) = kui*0.74_rkx*(log(z/zz0)- &
                      2.0_rkx*log((1.0_rkx+sqrt(1.0_rkx-9.0_rkx*zl))*0.5_rkx))
            end if
            ra(j,i) = max(ra(j,i),0.99_rkx)
            ra(j,i) = min(ra(j,i),999.9_rkx)
        end do
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine aerodyresis

    subroutine stomtresis(lai_f,laimin,laimax,ivegcov,igas, &
                          ustar,prec,sd,srad,ts,t2,rh,coszen,rc,rb)
      implicit none
      integer(ik4) , intent(in) :: igas
      integer(ik4) , intent(in) , dimension(jci1:jci2) :: ivegcov
      real(rkx) , dimension(jci1:jci2) , intent(in) :: coszen, srad , &
                         ts , rh , prec , sd , t2
      real(rkx) , dimension(jci1:jci2) , intent(in) :: lai_f , laimin , laimax
      real(rkx) , intent(in) , dimension(jci1:jci2) :: ustar
      real(rkx) , intent(out) , dimension(igas,jci1:jci2) :: rb , rc

      integer(ik4) :: j , l , lcov , ig
      real(rkx) :: rst, wst , rac , rgs_f
      real(rkx) :: rdu , rdv , rgo_f
      real(rkx) :: rcuto_f , rcuts_f
      real(rkx) :: ww1 , ww2 , ww3
      real(rkx) :: rdm , rdn , rv , rn
      real(rkx) :: ratio , sv , fv , fvv
      real(rkx) :: pardir , pardif
      real(rkx) :: tmaxk , tmink
      real(rkx) :: pshad , psun , rshad , rsun
      real(rkx) :: gshad , gsun , fsun , fshad
      real(rkx) :: gspar , temps !C
      real(rkx) :: bt , gt , gw , ryx
      real(rkx) :: es , d0 , gd , psi
      real(rkx) :: coedew , dq , usmin
      real(rkx) :: fsnow , rsnows
      real(rkx) :: dgas , di , vi
      real(rkx) :: dvh2o, rstom 
      real(rkx) :: rcut , rg , xp
      logical :: is_dew , is_rain
      real(rkx) , parameter :: dair = 0.369_rkx * 29.0_rkx + 6.29_rkx
      real(rkx) , parameter :: dh2o = 0.369_rkx * 18.0_rkx + 6.29_rkx
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'stomtresis'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif

      do j = jci1 , jci2
          is_rain = .false.
          is_dew  = .false.
          if ( ivegcov(j) == 0 ) then
            lcov = 14
          else
            lcov = ivegcov(j)
          end if
!         print*,'srad ====', srad(j)
!         print*,' ts  ====', ts(j)
!         print*,' coszen == ', coszen(j)

          tmaxk = tmax(lcov) + tzero
          tmink = tmin(lcov) + tzero
!         print *, ' tmax, tmin ==== ', tmaxk, tmink
          ! initialise rst as undef
          rst = -999.0_rkx
!         recalculate stomatal resistance according to Zhang et al.
          if (srad(j)   >= 0.1_rkx  .and. &
              ts(j)     <  tmaxk    .and. &
              ts(j)     >  tmink    .and. &
              lai_f(j)  > 0.001_rkx .and. &
              coszen(j) > 0.01_rkx ) then
            !================================================================
            ! Calculate direct and diffuse PAR from solar radiation and
            ! solar zenith angle
            !================================================================
            rdu   = 600.0_rkx * exp(-0.185_rkx/coszen(j))*coszen(j)
            rdv   = 0.4_rkx * (600.0_rkx - rdu ) * coszen(j)
            ww1   = -log(coszen(j))/2.302585_rkx
!           print *, ' ww1 = ', ww1
            ww2   = -1.195_rkx + 0.4459_rkx * ww1 - 0.0345_rkx * ww1**2
            ww3   = 1320.0_rkx*10.0_rkx**ww2
!           print *, 'ww= ', ww
            rdm   = (720.0_rkx*exp(-0.06_rkx/coszen(j))-ww3)*coszen(j)
!           print *, 'ww3= ', ww3, rdm
            rdn   = 0.6_rkx * (720.0_rkx - rdm - ww3) * coszen(j)
            rv    = max(0.1_rkx,  rdu + rdv)
            rn    = max(0.01_rkx, rdm + rdn)
            ratio = min(0.9_rkx,srad(j)/( rv + rn))
!           print *, 'ratio= ', ratio, rdn, rv, rn
            sv    = ratio * rv                            ! Total PAR
            fv    = min(0.99_rkx, (0.901_rkx - ratio)/0.7_rkx)
!           print *, 'sv  fv  = ', sv, fv
!           print *, 'rv  xxxxx  = ', rv, (1.0 - fv**0.6667)
            fvv   = max(0.01_rkx,rdu/rv*(1.0_rkx - fv**0.6667_rkx))
!           print *, 'fvv  = ', fvv
            ! fraction of PAR in the direct beam
            pardir = fvv * sv
            ! PAR from direct radiation
            pardif = sv - pardir
            ! PAR from diffuse radiation
!           print *, 'pardif========', pardif
            !================================================================
            ! Calculate sunlit and shaded leaf area, PAR for sunlit and
            ! shaded leaves
            !===============================================================
            if ( lai_f(j) > 2.5_rkx .and. srad(j) > 200.0_rkx ) then
              pshad = pardif * exp(-0.5_rkx * lai_f(j)**0.8_rkx) + &
                      0.07_rkx * pardir * (1.1_rkx-0.1_rkx*lai_f(j))* &
                      exp(-coszen(j))
              psun = pardir**0.8_rkx*0.5_rkx/coszen(j) + pshad
            else
              pshad = pardif * exp(-0.5_rkx * lai_f(j)**0.7_rkx) + &
                      0.07_rkx * pardir *(1.1_rkx-0.1_rkx*lai_f(j)) * &
                      exp(-coszen(j))
              psun = pardir * 0.5_rkx/coszen(j) + pshad
            end if
!           print *, 'pshad   psun   ', pshad , psun
            rshad = rsminz(lcov) + brs(lcov) * rsminz(lcov)/pshad
            rsun  = rsminz(lcov) + brs(lcov) * rsminz(lcov)/psun

            gshad = 1.0_rkx/rshad
            gsun  = 1.0_rkx/rsun
!           print *, 'rshad  ----< ', rshad, rsun, ' >---------rsun'
!           print *, 'gshad  ----< ', gshad, gsun, ' >-------- gsun'
            !================================================================
            ! Fsun, Fshade are the total sunlit and shaded leaf area
            ! index
            !================================================================
            xp = 0.5_rkx*lai_f(j)/coszen(j)
            if ( xp < 25.0_rkx ) then
              fsun  = 2.0_rkx*coszen(j)*(1.0_rkx-exp(-xp))
            else
              fsun = d_zero
            end if
            ! Sunlit leaf area
            fshad = lai_f(j) - fsun
            ! Shaded leaf area
!           print *, 'f, f ====',fshad,fsun
            !================================================================
            ! Stomatal conductance before including effects of
            ! temperature, vapor pressure defict and water stress.
            !================================================================
            gspar = fsun * gsun + fshad * gshad
            !================================================================
            ! function for temperature effect
            !================================================================
            temps = ts(j) - tzero
            bt = (tmax(lcov) - topt(lcov))/(topt(lcov) - tmin(lcov))
            gt = (tmax(lcov) - temps)/(tmax(lcov) - topt(lcov))
            gt = gt**bt
            gt = gt*(temps - tmin(lcov))/(topt(lcov) - tmin(lcov))
!           print *, 'gt ==========',gt
            !================================================================
            ! function for vapor pressure deficit
            !================================================================
            es = 6.108_rkx*exp(17.27_rkx*(ts(j)-tzero)/(ts(j)-35.86_rkx))
            d0 = es*(d_one-rh(j))/10.0_rkx ! kPa
            gd = 1.0_rkx - bvpd(lcov) * d0
!           print *, 'gd===',gd
            !================================================================
            ! function for water stress
            !================================================================
            psi = (-0.72_rkx - 0.0013_rkx * srad(j))
!           psi_s = (-0.395_rkx-0.043_rkx*(ts-tzero))*102.0_rkx
            gw = (psi - psi2(lcov))/(psi1(lcov) - psi2(lcov))
!           print *, 'gw==',gw
!           TEST
!           gw = 1
            if ( gw > 1.0_rkx ) gw = 1.0_rkx
            if ( gw < 0.1_rkx ) gw = 0.1_rkx
            if ( gd > 1.0_rkx ) gd = 1.0_rkx
            if ( gd < 0.1_rkx ) gd = 0.1_rkx
            !================================================================
            ! Stomatal resistance for water vapor
            !================================================================
            rst = 1.0_rkx / (gspar * gt * gd * gw)
!           print *, 'rst===',rst
          end if
          coedew = 0.1_rkx  ! for clear cloud
          es = 6.108_rkx*exp(17.27_rkx*(ts(j)-tzero)/(ts(j)-35.86_rkx))
          dq = 0.622_rkx/1000.0_rkx*es*(1.0_rkx-rh(j))*1000.0_rkx ! unit g/kg
          dq = max(0.0001_rkx,dq)
          usmin = 1.5_rkx/dq*coedew
!         print *, 'prec===== ', prec(j)
!         print *, 'usmin   ===  ', usmin
!         what is the unit of precipitation threshold
          if ( ts(j) > tzero .and. prec(j) > rainthr ) then
            is_rain = .true.
!           print *, 'rain==='
          else if (ts(j) > tzero .and. ustar(j) < usmin) then
            is_dew = .true.
!           print *, 'dew==='
!           print *, 'NO dew, NO rain ==='
          end if
!FAB TEST 
           !is_rain = .false.
           !is_dew = .false.

          !================================================================
          ! Decide fraction of stomatal blocking due to wet conditions
          !================================================================
          wst = 0.0_rkx
          if ( (is_dew .or. is_rain) .and. srad(j) > 200.0_rkx ) then
            wst = (srad(j) - 200.0_rkx)/800.0_rkx
            wst = min(wst, 0.5_rkx)
          end if
          !================================================================
          ! In-canopy aerodynamic resistance
          !================================================================
          rac = rac1(lcov)+(lai_f(j)-laimin(j))/ &
                (laimax(j)-laimin(j)+1.e-10_rkx)*(rac2(lcov)-rac1(lcov))
!         print *, 'rac1 = ', rac
          rac = rac*lai_f(j)**0.25_rkx/ustar(j)/ustar(j)
!         print *, 'rac2 = ', rac
          !================================================================
          ! Ground resistance for O3
          !================================================================
          if (ts(j) < 272.15_rkx .and. lcov /= 14 ) then
            rgo_f = min( rgo(lcov)*2.0_rkx, rgo(lcov) *     &
                           exp(0.2_rkx*(272.15_rkx-ts(j))))
!           print *, 'rgo_f1 =',rgo_f, ts(j)
          else
            rgo_f = rgo(lcov)
          end if
          !================================================================
          ! Ground resistance for SO2
          !================================================================
          if ( lcov == 12 ) then
            rgs_f = min(rgs(lcov)*(275.15_rkx - ts(j)), 500._rkx)
            rgs_f = max(rgs(lcov), 100._rkx)
!           print *, 'rgs_f ==== ', rgs_f
          else if ( is_rain .and. lcov /= 14 ) then
            rgs_f = 50.0_rkx
!           print *, 'rgs_f ==== ', rgs_f
          else if ( is_dew .and. lcov /= 14 ) then
            rgs_f = 100.0_rkx
!           print *, 'rgs_f ==== ', rgs_f
          else if ( ts(j) < 272.156_rkx .and. lcov /= 14 ) then
            rgs_f = min(rgs(lcov)*2.0_rkx, rgs(lcov) *     &
                          exp(0.2_rkx*(272.156_rkx - ts(j))))
!           print *, 'rgs_f ==== ', rgs_f
          else
            rgs_f = rgs(lcov)
!           print *, 'rgs_f ==== ', rgs_f
          end if
          !================================================================
          ! Cuticle resistance for O3 AND SO2
          !================================================================
          if ( rcutdo(lcov) <= -1.0_rkx ) then
            rcuto_f = 1.e25_rkx
            rcuts_f = 1.e25_rkx
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          else if ( is_rain ) then
            rcuto_f = rcutwo(lcov)/sqrt(lai_f(j))/ustar(j)
            rcuts_f = 50.0_rkx/sqrt(lai_f(j))/ustar(j)
            rcuts_f = max(rcuts_f, 20._rkx)
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          else if ( is_dew ) then
            rcuto_f = rcutwo(lcov)/sqrt(lai_f(j))/ustar(j)
            rcuts_f = 100.0_rkx/sqrt(lai_f(j))/ustar(j)
            rcuts_f = max(rcuts_f, 20._rkx)
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          else if (ts(j) < 272.156_rkx ) then
            ryx = exp(0.2_rkx * (272.156_rkx - ts(j) ))
            rcuto_f = rcutdo(lcov)/exp(3.0_rkx * rh(j))/     &
                      lai_f(j)**0.25_rkx/ustar(j)
            rcuts_f = rcutds(lcov)/exp(3.0_rkx * rh(j))/     &
                      lai_f(j)**0.25_rkx/ustar(j)
            rcuto_f = min(rcuto_f * 2.0_rkx, rcuto_f * ryx )
            rcuts_f = min(rcuts_f * 2.0_rkx, rcuts_f * ryx )
            rcuto_f = max(rcuto_f,100._rkx)
            rcuts_f = max(rcuts_f,100._rkx)
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          else
            rcuto_f = rcutdo(lcov)/exp(3.0_rkx*rh(j)) / &
                      lai_f(j)**0.25_rkx/ustar(j)
            rcuts_f = rcutds(lcov)/exp(3.0_rkx*rh(j)) / &
                      lai_f(j)**0.25_rkx/ustar(j)
            rcuto_f = max(rcuto_f, 100._rkx)
            rcuts_f = max(rcuts_f, 100._rkx)
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          end if
          !================================================================
          ! If snow occurs, Rg and Rcut are adjusted by snow cover
          ! fraction
          !================================================================
          fsnow = sd(j)/sdmax(lcov)
          fsnow = min(1.0_rkx, fsnow)   !snow cover fraction for leaves
!         print *, ' fsnow=  ', fsnow
          if ( fsnow > 0.0001_rkx .and. lcov /= 20 .or. &
                                        lcov /= 15 .or. &
                                        lcov /= 14 .or. &
                                        lcov /= 12 ) then
            rsnows = min(70.0_rkx*(275.15_rkx-ts(j)), 500._rkx)
            rsnows = max(rsnows, 100._rkx)
            rcuts_f = 1.0_rkx/((1.0_rkx - fsnow)/rcuts_f + fsnow/rsnows)
            rcuto_f = 1.0_rkx/((1.0_rkx - fsnow)/rcuto_f + fsnow/2000.0_rkx)
            fsnow = min(1.0_rkx, fsnow*2.0_rkx)
            ! snow cover fraction for ground
            rgs_f = 1.0_rkx/((1.0_rkx - fsnow)/rgs_f + fsnow/rsnows)
            rgo_f = 1.0_rkx/((1.0_rkx - fsnow)/rgo_f + fsnow/2000.0_rkx)
!           print *, 'rsnows= ', rsnows, ' fsnow=  ', fsnow
          end if
          !================================================================
          ! Calculate diffusivity for each gas species
          !================================================================
          do ig = 1 , igas
            dgas = 0.369_rkx * mw(ig) + 6.29_rkx
            di = 0.001_rkx*ts(j)**1.75_rkx * &
              sqrt((29.0_rkx + mw(ig))/mw(ig)/29._rkx)
            di = di/1.0_rkx/(dair**0.3333_rkx + dgas**0.3333_rkx)**2
            vi = 145.8_rkx * 1.e-4_rkx * &
                 (ts(j) * 0.5_rkx + t2(j) *0.5_rkx)**1.5_rkx/ &
                 (ts(j) * 0.5_rkx + t2(j) *0.5_rkx + 110.4_rkx)
            !================================================================
            ! Calculate quasi-laminar resistance
            !================================================================
            rb(ig,j) = 5.0_rkx/ustar(j) * (vi/di)**.666667_rkx
!           print *, 'rb==', rb(ig,l,j)
            !================================================================
            ! Calculate stomatal resistance for each species from the ratio
            ! of  diffusity of water vapor to the gas species
            !================================================================
            dvh2o = 0.001_rkx*ts(j)**1.75_rkx * &
              sqrt((29.0_rkx+18.0_rkx)/29.0_rkx/18.0_rkx)
            dvh2o = dvh2o/(dair**0.3333_rkx + dh2o**0.3333_rkx)**2
            rstom = rst * dVh2o/di + rm(ig)
            ! (rst <999) for bare surfaces)
            !================================================================
            ! Scale cuticle and ground resistances for each species
            !================================================================
            rcut = 1.0_rkx/(alphaz(ig)/rcuts_f+betaz(ig)/rcuto_f)
            rg   = 1.0_rkx/(alphaz(ig)/rgs_f+betaz(ig)/rgo_f)
            !================================================================
            ! Calculate total surface resistance
            !================================================================
            ! account for zero stomatal resistance (rst and rstom are zero
            ! for bare surfaces)
            ! set wst to 1 also in that case (total stomatal blocking).
            if ( rst == -999.0 ) wst = 1.0_rkx
!           rc(ig,l,j) = (1.0_rkx - wst)/rstom + 1.0_rkx/(rg)+1.0_rkx/rcut
            rc(ig,j) = (1.0_rkx - wst)/rstom + 1.0_rkx/(rac+rg)+1.0_rkx/rcut
            rc(ig,j) = max(10._rkx,1.0_rkx/rc(ig,j))
          end do !igas
      end do !ilg
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine stomtresis

end module mod_che_drydep
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
