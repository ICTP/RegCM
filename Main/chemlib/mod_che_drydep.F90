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
  real(rkx) , parameter :: a1 = 145.8_rkx
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
  ! Only one cover type per grid cell for now
  !
  integer, parameter :: luc = 1
  !
  ! Number of gas taken into account by drydep scheme
  !
  integer, parameter :: ngasd = 31
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

    subroutine drydep_aero(j,mbin,indsp,rhop,ivegcov,throw,roarow, &
                           shj,pressg,temp2,sutemp,srad,rh10,      &
                           wind10,zeff,beffdiam,pdepv,ddepv)
      implicit none
      integer(ik4) , intent(in) :: j , mbin
      integer(ik4) , intent(in) , dimension(mbin) :: indsp
      integer(ik4) , intent(in) , dimension(ici1:ici2) :: ivegcov
      real(rkx) , dimension(ici1:ici2) , intent(in) :: pressg , rh10 , &
                       srad , sutemp , temp2 , wind10 , zeff
      real(rkx) , dimension(ici1:ici2,kz) , intent(in) :: roarow , throw
      real(rkx) , dimension(kz) , intent(in) :: shj
      real(rkx) , dimension(mbin) , intent(in) :: beffdiam
      real(rkx) , intent(in) :: rhop

      ! output table to be passed out. Care dimension is ntr

      real(rkx) , intent(out) , dimension(ici1:ici2,kz,ntr) :: pdepv
      real(rkx) , intent(out) , dimension(ici1:ici2,ntr) :: ddepv

      real(rkx) :: amfp , amob , eb , eim , ein , frx1
      real(rkx) :: pre , prii , priiv , r1 , st
      real(rkx) , dimension(ici1:ici2,kz) :: amu
      real(rkx) , dimension(ici1:ici2) :: anu , schm
      real(rkx) , dimension(ici1:ici2,kz,mbin) :: cfac , pdepvsub , pdiff , &
                  rhsize , taurel
      real(rkx) , dimension(ici1:ici2,luc) :: ra , ustar
      real(rkx) , dimension(ici1:ici2,luc,mbin) :: rs
      real(rkx), dimension(ici1:ici2,2:kz) :: wk, settend
      real(rkx) , dimension(mbin) :: avesize
      integer(ik4) :: i , k , kcov , l , n , ib
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'drydep_aero'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      ! here avesize is a RADIUS of the particle bin in m :
      ! calculated here from bin effective diameter in micrometer
      do n = 1 , mbin
        avesize(n) = beffdiam(n)* 1.e-6_rkx * d_half
      end do
      ! ********************************************************
      ! *   aerosize - dry radius                    !!     ****
      ! *   rhop  - density for each aerosol type           ****
      ! ********************************************************
      do n = 1 , mbin
        do l = 1 , kz
          do i = ici1 , ici2
            !
            ! ********************************************************
            ! *  aerosol gravitational settling velocity          ****
            ! *  and diffusion coefficient                        ****
            ! *                                                   ****
            ! * air's dynamic viscosity                           ****
            ! ********************************************************
            !
            amu(i,l) = a1*1.e-8_rkx*throw(i,l)**a2/(throw(i,l)+a3)
            ! mid layer pressure in [pascal].
            pre = pressg(i)*shj(l)
            !
            ! ********************************************************
            ! * mean molecular free path.                         ****
            ! *     k.v. beard [1976], j atm. sci., 33            ****
            ! ********************************************************
            !
            amfp = c1*(amu(i,l)/c2)*(c3/pre)*sqrt(throw(i,l)/c4)
            prii = 2.0_rkx/9.0_rkx*egrav/amu(i,l)
            priiv = prii*(rhop-roarow(i,l))
            !
            ! ********************************************************
            ! * cunningham slip correction factor and             ****
            ! * relaxation time = vg/grav.                        ****
            ! ********************************************************
            !
            cfac(i,l,n) = 1.0_rkx + amfp/avesize(n) * &
                         (aa1+aa2*exp(-aa3*avesize(n)/amfp))
            taurel(i,l,n) = max(priiv*avesize(n)**2*cfac(i,l,n) * &
                            regrav,0._rkx)
            !
            ! ********************************************************
            ! * stokes friction                                  *****
            ! ! pdepvsub(i,l,n) ' sellting dep. velocity = '
            ! ********************************************************
            !
            pdepvsub(i,l,n) = taurel(i,l,n)*egrav
          end do
        end do
      end do
      !
      ! Find aerodynamic resistance
      !
      call aerodyresis(zeff,wind10,temp2,sutemp,rh10,srad,ivegcov,ustar,ra)
      !
      ! *****************************************************
      ! * the schmidt number is the ratio of the         ****
      ! * kinematic viscosity of air to the particle     ****
      ! * brownian diffusivity ===> sc=v/d               ****
      ! *****************************************************
      !
      do n = 1 , mbin
        do l = 1 , kz
          do i = ici1 , ici2
            !
            ! *****************************************************
            ! * for now we will not consider the humidity      ****
            ! * impact so we will set the variable frx1=1.0    ****
            ! * i.e. only dry particles                        ****
            ! *****************************************************
            !
            frx1 = 1.0_rkx
            rhsize(i,l,n) = avesize(n)*frx1 ! still a radius
            anu(i) = amu(i,l)/roarow(i,l)
            amob = 6.0_rkx*mathpi*amu(i,l)*rhsize(i,l,n)/cfac(i,l,n)
            pdiff(i,l,n) = boltzk*throw(i,l)/amob
            schm(i) = anu(i)/pdiff(i,l,n)
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
          if ( l == kz ) then
            do k = 1 , luc ! luc  = 1 for the moment
              do i = ici1 , ici2
                !
                ! find the right table index for the cell cover ( ocean
                ! and lake are 0 in the ivegcov and 14-15 in the table )
                !
                if ( ivegcov(i) == 0 ) then
                  kcov = 14
                else if ( ivegcov(i) > 20 ) then
                  kcov = 20
                else
                  kcov = ivegcov(i)
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
                st = taurel(i,l,n)*ustar(i,k)*ustar(i,k)/anu(i)
                eb = schm(i)**(-0.666667_rkx)
!               eim=(st/(st+aest(k)))**2
                eim = (st/(st+aest(kcov)))**2
                eim = min(eim,0.6_rkx)
                ein = 0.0_rkx
!               if (arye(k) > 0.0001_rkx) then
!                 ein = (1000.0_rkx*2.0_rkx*avesize(n)/arye(k))**1.5_rkx
!               end if
                if ( arye(kcov) > 0.0001_rkx ) then
                  ein = (1000.0_rkx*2.0_rkx*avesize(n)/arye(kcov))**1.5_rkx
                end if
                ein = min(ein,0.5_rkx)
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

                r1 = max(0.5_rkx,exp(-(min(sqrt(st),25.0_rkx))))
                if ( kcov >= 11 .and. r1 < 0.5_rkx ) r1 = 0.5_rkx
                if ( r1 < 0.4_rkx ) r1 = 0.4_rkx
                ! ***************************************************
                ! * calculation of rs: the surface resistance   *****
                ! * which depends on the collection efficiency  *****
                ! * of the surface and is determined by the     *****
                ! * various deposition processes                *****
                ! ***************************************************
!               rs= 1.0/ustar(i,k)/(eb+eim+ein)/r1
                rs(i,k,n) = 1.0_rkx/3.0_rkx/ustar(i,k)/(eb+eim+ein)/r1
              end do
            end do
          end if
        end do
      end do

      ! average settling and deposition velocities on bin
      ! care we use pdepv and ddpv table that are dimensionned to ntr
      ! and not mbin !
      do ib = 1 , mbin
        ! there isw no sub-bin anymore / we consider directly effective radius
        pdepv(:,:,indsp(ib)) = 0.0_rkx
        ddepv(:,indsp(ib))   = 0.0_rkx
        do i = ici1 , ici2
          pdepv(i,:,indsp(ib)) = pdepvsub(i,:,ib)
          ! agregate the dry deposition velocity, remember one cover per grid
          ! cell for now
          ! the dry deposition deposition velocity must accound also for the
          ! settling vrlocity at kz
          ! simple form now add the vs
          ddepv(i,indsp(ib)) = 1.0_rkx/(ra(i,1)+rs(i,1,ib)) + pdepvsub(i,kz,ib)
        end do
      end do
      !
      ! Finally update the emission and settling tendencies for
      ! dust and sea salt
      !
      do ib = 1 , mbin
        ! deposition, remember chiten must be normalised by psb and
        ! consistent with chib
        do k = 2 , kz
          do i = ici1 , ici2
            wk(i,k) =  max(twt(k,1)*chib(j,i,k,indsp(ib)) + &
                           twt(k,2)*chib(j,i,k-1,indsp(ib)),mintr)
          end do
        end do
        do i = ici1 , ici2
          do k = 2 , kz - 1
            ! do not apply to the first level
            !settend(i,k) = (wk(i,k+1)*pdepv(i,k+1,indsp(ib)) - &
            !                wk(i,k)*pdepv(i,k,indsp(ib))) / cdzq(j,i,k)
            ! use exponential form for stability
            settend(i,k) =  wk(i,k+1) * (d_one - &
              exp(-ddepv(i,indsp(ib))/cdzq(j,i,k)*dt ))/dt  - &
                            wk(i,k)   * (d_one - &
              exp(-ddepv(i,indsp(ib))/cdzq(j,i,k)*dt ))/dt

            chiten(j,i,k,indsp(ib)) = chiten(j,i,k,indsp(ib)) - settend(i,k)
            if ( ichdiag == 1 ) then
              cseddpdiag(j,i,k,indsp(ib)) = cseddpdiag(j,i,k,indsp(ib)) - &
                                            settend(i,k) * cfdout
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

            ! settend(i,kz) =  (chib(j,i,kz,indsp(ib))  * ddepv(i,indsp(ib))-  &
            !                   wk(i,kz)*pdepv(i,kz,indsp(ib))) / cdzq(j,i,kz)
            settend(i,kz) = max(chib(j,i,kz,indsp(ib)),mintr) * &
              (d_one - exp(-ddepv(i,indsp(ib))/cdzq(j,i,kz)*dt ))/dt -  &
                            wk(i,kz) * &
              (d_one - exp(-pdepv(i,kz,indsp(ib))/cdzq(j,i,k)*dt))/dt
            chiten(j,i,kz,indsp(ib)) = chiten(j,i,kz,indsp(ib)) - settend(i,kz)

            ! save the dry deposition flux for coupling with
            ! landsurface scheme (Kg.m2.s-1)
            ! consider ddflux = Cav . Vd where Cav would be the average
            ! concentration within the time step Cav = 0.5 (C + (C+deltaC))
            ! care  chib and settend have to be corrected for pressure
            ! cdrydepflux is a time accumulated array set to zero when surface
            ! scheme is called (cf atm to surf interface)
            cdrydepflx(j,i,indsp(ib)) = cdrydepflx(j,i,indsp(ib)) + &
              (chib(j,i,kz,indsp(ib))  - settend(i,kz)*dt/d_two) / &
                 cpsb(j,i) * crhob3d(j,i,kz) * ddepv(i,indsp(ib))

            !diagnostic for settling and drydeposition removal
            if ( ichdiag == 1 ) then
              cseddpdiag(j,i,kz,indsp(ib)) = cseddpdiag(j,i,kz,indsp(ib)) - &
                                             settend(i,kz) * cfdout
            end if

            ! accumulated diagnostic for dry deposition flux
            ! average (in kg .m2.s-1)
            remdrd(j,i,indsp(ib)) = remdrd(j,i,indsp(ib)) + &
                                    cdrydepflx(j,i,indsp(ib)) * cfdout

            ! no net flux is passed to BL schemes in this case
            chifxuw(j,i,indsp(ib)) = d_zero
            drydepv(j,i,indsp(ib)) = d_zero

          else if ( ichdrdepo == 2 ) then
            !
            ! add the dry deposition term to the net emision/deposition flux
            ! for the BL scheme !
            ! flux
            chifxuw(j,i,indsp(ib)) = chifxuw(j,i,indsp(ib)) - &
                chib(j,i,kz,indsp(ib))/ cpsb(j,i) * ddepv(i,indsp(ib))

            drydepv(j,i,indsp(ib)) = ddepv(i,indsp(ib))

          end if
          !
          ! dry dep velocity diagnostic in m.s-1  ( + drydep v. include
          ! also settling , accumulated between two outputs time step)
          ddv_out(j,i,indsp(ib)) = ddv_out(j,i,indsp(ib)) + &
                 ddepv(i,indsp(ib))
        end do
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine drydep_aero

    subroutine drydep_gas(j,lmonth,lday,ivegcov,rh10,srad,tsurf, &
                          prec,temp10,wind10,zeff)
      use mod_che_indices
      implicit none
      integer(ik4) , intent(in) :: j
      integer, intent(in) :: lmonth , lday
      integer(ik4) , intent(in) , dimension(ici1:ici2) :: ivegcov
      real(rkx) , intent(in) , dimension(ici1:ici2) :: rh10 , srad , tsurf , &
                                            prec, temp10 , wind10 , zeff
      real(rkx),  dimension(ici1:ici2,ntr) :: drydepvg

      integer(ik4) :: n , i , im , kcov
      real(rkx) , dimension(ici1:ici2,luc) :: ustar, resa
      real(rkx) , dimension(ngasd,ici1:ici2,luc) :: resb, resc
      real(rkx) , dimension(ngasd,ici1:ici2,luc) :: vdg
      real(rkx) , dimension(ici1:ici2) :: icz , ddrem
      real(rkx) , dimension(ici1:ici2) :: lai_f , laimin , laimax , snow
      real(rkx) :: kd , kav , rdz
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'drydep_gas'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      ! Different options for LAI and roughness
      ! for the moment read from

      do i = ici1 , ici2
        if ( ivegcov(i) == 0 ) then
          kcov = 14
        else if ( ivegcov(i) > 20 ) then
          kcov = 20
        else
          kcov = ivegcov(i)
        end if
        im = lmonth - 1
        if ( lmonth == 1 ) im = 12
        if (lday <= 15 ) then
          lai_f(i) = lai(kcov,im) + (lai(kcov,lmonth) - &
                     lai(kcov,im))/30._rkx * real(15 + lday,rkx)
        else
          lai_f(i) = lai(kcov,lmonth) + (lai(kcov,lmonth+1) - &
                     lai(kcov,lmonth))/30._rkx * real(lday - 15,rkx)
        end if
        if ( lai_f(i) < d_zero) lai_f(i) = d_zero
        laimin(i) = lai(kcov,14)
        laimax(i) = lai(kcov,15)
      end do
      call aerodyresis(zeff,wind10,temp10,tsurf,rh10,srad,ivegcov,ustar,resa)
      snow(:) = d_zero
      icz(:) = czen(j,:)
      call stomtresis(lai_f,laimin,laimax,ivegcov,ngasd,ustar,prec,snow,srad, &
                      tsurf,temp10,rh10,icz,resc,resb)
      ! now calculate the dry deposition velocities and select it
      ! according to the gasphase mechanism
      ! vdg in m.s-1
       vdg(:,:,:) = d_zero
       do i = ici1 , ici2
         do n = 1 , ngasd
           vdg(n,i,:) = d_one/(resa(i,:)+resb(n,i,:)+resc(n,i,:))
         end do
       end do
       ! this part depends on the chem mechanism
       ! for CBMZ , we can certainly improve this.

       drydepvg = d_zero
       drydepvg(:,iso2)  =  vdg(1,:,1)
       ! SO2 deposition is used in SULF , AERO and CBMZ simulations
       if ( igaschem > 0 ) then
         drydepvg(:,ino2)  =  vdg(3,:,1)!*0.5
         drydepvg(:,io3)   =  vdg(4,:,1)!*0.5
         drydepvg(:,ih2o2) =  vdg(5,:,1)!*0.5
         drydepvg(:,ihno3) =  vdg(6,:,1)!*0.5
!        drydepvg(:,inh3)  =  vdg(9,:,1)!*0.5
         drydepvg(:,ipan)  =  vdg(10,:,1)!*0.5
         drydepvg(:,ihcho) =  vdg(14,:,1)!*0.5
         drydepvg(:,iald2) =  vdg(15,:,1)!*0.5
         drydepvg(:,ich3oh)  =  vdg(23,:,1)!*0.5
       end if

       ! Finally : gas phase dry dep tendency calculation
       if ( ichdrdepo == 1 ) then
         do i = ici1 , ici2
           rdz = d_one / cdzq(j,i,kz)

           ! If using CLM then use the dry deposition velocities coming directly
           ! from internal CLM calculations
#ifdef CLM
           if ( ivegcov(i) == 0 ) then
             drydepvg(i,:) = cdep_vels(j,i,:)
           end if
#endif
           do n = 1 , ntr
             kd = drydepvg(i,n) * rdz !Kd removal rate in s-1
             kav = max(chib(j,i,kz,n)-mintr,d_zero)
             if ( kd*dt < 25.0_rkx ) then
               ! dry dep removal tendency (+)
               ddrem(i) = kav * (d_one-exp(-kd*dt))/dt
             else
               ddrem(i) = d_zero
             end if
             ! update chiten
             chiten(j,i,kz,n) = chiten(j,i,kz,n) - ddrem(i)
             ! diag dry dep tendency
             if ( ichdiag == 1 ) then
                cseddpdiag(j,i,kz,n) = cseddpdiag(j,i,kz,n) - &
                                                ddrem(i) * cfdout
             end if
             ! drydep flux diagnostic (accumulated between two outputs time
             ! step) ! flux is in kg/m2/s-1 so need to normalise by ps here.
             remdrd(j,i,n) = remdrd(j,i,n) + ddrem(i)/cpsb(j,i) * cfdout
             ! dry dep velocity diagnostic in m.s-1
             ! (accumulated between two outputs time step)
             drydepv(j,i,n) = d_zero
             ddv_out(j,i,n) = ddv_out(j,i,n) + drydepvg(i,n)
           end do
         end do
       else if ( ichdrdepo == 2 ) then
         do i = ici1 , ici2
           ! If using CLM then use the dry deposition velocities coming directly
           ! from internal CLM calculations
#ifdef CLM
           if ( ivegcov(i) == 0 ) then
             drydepvg(i,:)  =  cdep_vels(j,i,:)
           end if
#endif
           do n = 1 , ntr
             chifxuw(j,i,n) = chifxuw(j,i,n) - (chib(j,i,kz,n) / &
                                 cpsb(j,i)) * drydepvg(i,n)
             ! dry dep velocity diagnostic in m.s-1
             ! (accumulated between two outputs time step)
             drydepv(j,i,n) =  drydepvg(i,n)
             ddv_out(j,i,n) =  ddv_out(j,i,n) + drydepvg(i,n)
           end do
         end do

       end if
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine drydep_gas

    subroutine aerodyresis(zeff,wind10,temp2,sutemp,rh10,srad,ivegcov,ustar,ra)
      implicit none
      integer(ik4) , dimension(ici1:ici2) , intent(in) :: ivegcov
      real(rkx) , dimension(ici1:ici2) , intent(in) :: temp2 , wind10 , rh10
      real(rkx) , dimension(ici1:ici2) , intent(in) :: sutemp , srad , zeff
      real(rkx) , dimension(ici1:ici2,luc) , intent(out) :: ustar , ra
      integer(ik4) :: i , j
      real(rkx) :: vp , tsv
      real(rkx) :: z , zl , ww
      real(rkx) :: ptemp2 , es , qs
      real(rkx) :: wvpm , vptemp , tsw , mol
      real(rkx) :: z0water , dthv , cun , zdl
      real(rkx) :: psiu , psit , x , y
      real(rkx) :: thstar , rib , dtemp , tbar
      real(rkx) :: ustarsq , utstar , kui
      real(rkx) :: ratioz , logratio , asq
      real(rkx) :: aa , cm , ch , fm , fh
      real(rkx) , dimension(ici1:ici2) :: zz0
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
      do j = 1 , luc
        do i = ici1 , ici2
          ww = max(wind10(i),1.0_rkx)
          zz0(i) = zeff(i)
          ! ***************************************************************
          ! * potential temperature at z2  (deg. k)
          ! ***************************************************************
          ptemp2 = temp2(i) + z10*0.0098_rkx
          ! ***************************************************************
          ! * for calculations over water compute values of critical
          ! * profile variables: l and ustar
          ! *           ******begin for water***
          ! ***************************************************************
          if ( ivegcov(i) == 0 ) then
            ! **************************************************************
            ! * vp  - vapour pressure at z2
            ! * wvpm- water vapour mixing ratio at  z2
            ! * vptemp- virtual potential temperature at z2
            ! **************************************************************
            es = 6.108_rkx*exp(17.27_rkx*(temp2(i)-tzero)/(temp2(i)-35.86_rkx))
            vp = rh10(i)*es
            wvpm = ep2*vp/(stdpmb-vp)
            vptemp = ptemp2*(1.0_rkx+0.61_rkx*wvpm)
            ! **************************************************************
            ! * assume rh10 at water surface is 100%
            ! *   vp = es(tsw-tzero) !sat. vap press at surface
            ! *   saturated vapour pressure at surface
            ! *   saturated mixing ratio at surface
            ! *   tsv - virtual potential temperature at surface (deg. k)
            ! **************************************************************
            tsw = sutemp(i)
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
              psiu = 2.0_rkx*log(0.5_rkx*(1.0_rkx+x)) + \
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
            ustar(i,j) = vonkar*ww/(log(z10/z0water)-psiu)
            thstar = vonkar*(ptemp2-sutemp(i)) / &
                     (0.74_rkx*log(z10/z0water)-psit)
            zz0(i) = z0water
          else
            ! **************************************************************
            ! * compute ustar and l for land use categories other than
            ! * water use louis method. !pkk 7/16/85, find bulk
            ! * richardson number.
            ! **************************************************************
            rib = egrav*z10*(ptemp2-sutemp(i))/(sutemp(i)*ww**2)
            ! ***************************************************************
            ! * ensure that conditions over land are never stable when
            ! * there is incoming solar radiation
            ! ***************************************************************
            if ( srad(i) > 0.0_rkx .and. rib > 0.0_rkx ) rib = 1.e-15_rkx
            dtemp = ptemp2 - sutemp(i)
            if ( abs(dtemp) < 1.e-10_rkx ) dtemp = sign(1.e-10_rkx,dtemp)
            tbar = 0.5_rkx*(ptemp2+sutemp(i))
            ratioz = z10/zz0(i)
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
            ustar(i,j) = sqrt(ustarsq)
            thstar = utstar/ustar(i,j)
            mol = tbar*ustarsq/(vonkar*egrav*thstar)
          end if

          kui = 1.0_rkx/(vonkar*ustar(i,j))

          ! **************************************************************
          ! * compute the values of  ra                            *******
          ! **************************************************************
          z = z10
          zl = z/mol
          if ( zl >= 0.0_rkx ) then
            ra(i,j) = kui*(0.74_rkx*log(z/zz0(i))+4.7_rkx*zl)
          else
            ra(i,j) = kui*0.74_rkx*(log(z/zz0(i))- &
                      2.0_rkx*log((1.0_rkx+sqrt(1.0_rkx-9.0_rkx*zl))*0.5_rkx))
          end if
          ra(i,j) = max(ra(i,j),0.99_rkx)
          ra(i,j) = min(ra(i,j),999.9_rkx)
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
      integer(ik4) , intent(in) , dimension(ici1:ici2) :: ivegcov
      real(rkx) , dimension(ici1:ici2) , intent(in) :: coszen, srad , &
                         ts , rh , prec , sd , t2
      real(rkx) , dimension(ici1:ici2) , intent(in) :: lai_f , laimin , laimax
      real(rkx) , intent(in) , dimension(ici1:ici2,luc) :: ustar
      real(rkx) , intent(out) , dimension(igas,ici1:ici2,luc) :: rb , rc

      integer(ik4) :: i , j , kcov , ig
      real(rkx) :: rst , wst , rac , rgs_f
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
      real(rkx) :: dvh2o , rstom
      real(rkx) :: rcut , rg , xp
      logical :: is_dew , is_rain
      real(rkx) , parameter :: dair = 0.369_rkx * 29.0_rkx + 6.29_rkx
      real(rkx) , parameter :: dh2o = 0.369_rkx * 18.0_rkx + 6.29_rkx
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'stomtresis'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif

      do j = 1 , luc
        do i = ici1 , ici2
          is_rain = .false.
          is_dew  = .false.
          if ( ivegcov(i) == 0 ) then
            kcov = 14
          else
            kcov = ivegcov(i)
          end if
!         print*,'srad ====', srad(i)
!         print*,' ts  ====', ts(i)
!         print*,' coszen == ', coszen(i)

          tmaxk = tmax(kcov) + tzero
          tmink = tmin(kcov) + tzero
!         print *, ' tmax, tmin ==== ', tmaxk, tmink
          ! initialise rst as undef
          rst = -999.0_rkx
          if (srad(i)   >= 0.1_rkx  .and. &
              ts(i)     <  tmaxk    .and. &
              ts(i)     >  tmink    .and. &
              lai_f(i)  > 0.001_rkx .and. &
              coszen(i) > 0.01_rkx ) then
            !================================================================
            ! Calculate direct and diffuse PAR from solar radiation and
            ! solar zenith angle
            !================================================================
            rdu   = 600.0_rkx * exp(-0.185_rkx/coszen(i))*coszen(i)
            rdv   = 0.4_rkx * (600.0_rkx - rdu ) * coszen(i)
            ww1   = -log(coszen(i))/2.302585_rkx
!           print *, ' ww1 = ', ww1
            ww2   = -1.195_rkx + 0.4459_rkx * ww1 - 0.0345_rkx * ww1**2
            ww3   = 1320.0_rkx*10.0_rkx**ww2
!           print *, 'ww= ', ww
            rdm   = (720.0_rkx*exp(-0.06_rkx/coszen(i))-ww3)*coszen(i)
!           print *, 'ww3= ', ww3, rdm
            rdn   = 0.6_rkx * (720.0_rkx - rdm - ww3) * coszen(i)
            rv    = max(0.1_rkx,  rdu + rdv)
            rn    = max(0.01_rkx, rdm + rdn)
            ratio = min(0.9_rkx,srad(i)/( rv + rn))
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
            if ( lai_f(i) > 2.5_rkx .and. srad(i) > 200.0_rkx ) then
              pshad = pardif * exp(-0.5_rkx * lai_f(i)**0.8_rkx) + &
                      0.07_rkx * pardir * (1.1_rkx-0.1_rkx*lai_f(i))* &
                      exp(-coszen(i))
              psun = pardir**0.8_rkx*0.5_rkx/coszen(i) + pshad
            else
              pshad = pardif * exp(-0.5_rkx * lai_f(i)**0.7_rkx) + &
                      0.07_rkx * pardir *(1.1_rkx-0.1_rkx*lai_f(i)) * &
                      exp(-coszen(i))
              psun = pardir * 0.5_rkx/coszen(i) + pshad
            end if
!           print *, 'pshad   psun   ', pshad , psun
            rshad = rsminz(kcov) + brs(kcov) * rsminz(kcov)/pshad
            rsun  = rsminz(kcov) + brs(kcov) * rsminz(kcov)/psun
            gshad = 1.0_rkx/rshad
            gsun  = 1.0_rkx/rsun
!           print *, 'rshad  ----< ', rshad, rsun, ' >---------rsun'
!           print *, 'gshad  ----< ', gshad, gsun, ' >-------- gsun'
            !================================================================
            ! Fsun, Fshade are the total sunlit and shaded leaf area
            ! index
            !================================================================
            xp = 0.5_rkx*lai_f(i)/coszen(i)
            if ( xp < 25.0_rkx ) then
              fsun  = 2.0_rkx*coszen(i)*(1.0_rkx-exp(-xp))
            else
              fsun = d_zero
            end if
            ! Sunlit leaf area
            fshad = lai_f(i) - fsun
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
            temps = ts(i) - tzero
            bt = (tmax(kcov) - topt(kcov))/(topt(kcov) - tmin(kcov))
            gt = (tmax(kcov) - temps)/(tmax(kcov) - topt(kcov))
            gt = gt**bt
            gt = gt*(temps - tmin(kcov))/(topt(kcov) - tmin(kcov))
!           print *, 'gt ==========',gt
            !================================================================
            ! function for vapor pressure deficit
            !================================================================
            es = 6.108_rkx*exp(17.27_rkx*(ts(i)-tzero)/(ts(i)-35.86_rkx))
            d0 = es*(d_one-rh(i))/10.0_rkx ! kPa
            gd = 1.0_rkx - bvpd(kcov) * d0
!           print *, 'gd===',gd
            !================================================================
            ! function for water stress
            !================================================================
            psi = (-0.72_rkx - 0.0013_rkx * srad(i))
!           psi_s = (-0.395_rkx-0.043_rkx*(ts-tzero))*102.0_rkx
            gw = (psi - psi2(kcov))/(psi1(kcov) - psi2(kcov))
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
          es = 6.108_rkx*exp(17.27_rkx*(ts(i)-tzero)/(ts(i)-35.86_rkx))
          dq = 0.622_rkx/1000.0_rkx*es*(1.0_rkx-rh(i))*1000.0_rkx ! unit g/kg
          dq = max(0.0001_rkx,dq)
          usmin = 1.5_rkx/dq*coedew
!         print *, 'prec===== ', prec(i)
!         print *, 'usmin   ===  ', usmin
!         what is the unit of precipitation threshold
          if ( ts(i) > tzero .and. prec(i) > rainthr ) then
            is_rain = .true.
!           print *, 'rain==='
          else if (ts(i) > tzero .and. ustar(i,j) < usmin) then
            is_dew = .true.
!           print *, 'dew==='
!           print *, 'NO dew, NO rain ==='
          end if
          !================================================================
          ! Decide fraction of stomatal blocking due to wet conditions
          !================================================================
          wst = 0.0_rkx
          if ( (is_dew .or. is_rain) .and. srad(i) > 200.0_rkx ) then
            wst = (srad(i) - 200.0_rkx)/800.0_rkx
            wst = min(wst, 0.5_rkx)
          end if
          !================================================================
          ! In-canopy aerodynamic resistance
          !================================================================
          rac = rac1(kcov)+(lai_f(i)-laimin(i))/ &
                (laimax(i)-laimin(i)+1.e-10_rkx)*(rac2(kcov)-rac1(kcov))
!         print *, 'rac1 = ', rac
          rac = rac*lai_f(i)**0.25_rkx/ustar(i,j)/ustar(i,j)
!         print *, 'rac2 = ', rac
          !================================================================
          ! Ground resistance for O3
          !================================================================
          if (ts(i) < 272.15_rkx .and. kcov /= 14 ) then
            rgo_f = min( rgo(kcov)*2.0_rkx, rgo(kcov) *     &
                           exp(0.2_rkx*(272.15_rkx-ts(i))))
!           print *, 'rgo_f1 =',rgo_f, ts(i)
          else
            rgo_f = rgo(kcov)
          end if
          !================================================================
          ! Ground resistance for SO2
          !================================================================
          if ( kcov == 12 ) then
            rgs_f = min(rgs(kcov)*(275.15_rkx - ts(i)), 500._rkx)
            rgs_f = max(rgs(kcov), 100._rkx)
!           print *, 'rgs_f ==== ', rgs_f
          else if ( is_rain .and. kcov /= 14 ) then
            rgs_f = 50.0_rkx
!           print *, 'rgs_f ==== ', rgs_f
          else if ( is_dew .and. kcov /= 14 ) then
            rgs_f = 100.0_rkx
!           print *, 'rgs_f ==== ', rgs_f
          else if ( ts(i) < 272.156_rkx .and. kcov /= 14 ) then
            rgs_f = min(rgs(kcov)*2.0_rkx, rgs(kcov) *     &
                          exp(0.2_rkx*(272.156_rkx - ts(i))))
!           print *, 'rgs_f ==== ', rgs_f
          else
            rgs_f = rgs(kcov)
!           print *, 'rgs_f ==== ', rgs_f
          end if
          !================================================================
          ! Cuticle resistance for O3 AND SO2
          !================================================================
          if ( rcutdo(kcov) <= -1.0_rkx ) then
            rcuto_f = 1.e25_rkx
            rcuts_f = 1.e25_rkx
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          else if ( is_rain ) then
            rcuto_f = rcutwo(kcov)/sqrt(lai_f(i))/ustar(i,j)
            rcuts_f = 50.0_rkx/sqrt(lai_f(i))/ustar(i,j)
            rcuts_f = max(rcuts_f, 20._rkx)
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          else if ( is_dew ) then
            rcuto_f = rcutwo(kcov)/sqrt(lai_f(i))/ustar(i,j)
            rcuts_f = 100.0_rkx/sqrt(lai_f(i))/ustar(i,j)
            rcuts_f = max(rcuts_f, 20._rkx)
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          else if (ts(i) < 272.156_rkx ) then
            ryx = exp(0.2_rkx * (272.156_rkx - ts(i) ))
            rcuto_f = rcutdo(kcov)/exp(3.0_rkx * rh(i))/     &
                      lai_f(i)**0.25_rkx/ustar(i,j)
            rcuts_f = rcutds(kcov)/exp(3.0_rkx * rh(i))/     &
                      lai_f(i)**0.25_rkx/ustar(i,j)
            rcuto_f = min(rcuto_f * 2.0_rkx, rcuto_f * ryx )
            rcuts_f = min(rcuts_f * 2.0_rkx, rcuts_f * ryx )
            rcuto_f = max(rcuto_f,100._rkx)
            rcuts_f = max(rcuts_f,100._rkx)
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          else
            rcuto_f = rcutdo(kcov)/exp(3.0_rkx*rh(i)) / &
                      lai_f(i)**0.25_rkx/ustar(i,j)
            rcuts_f = rcutds(kcov)/exp(3.0_rkx*rh(i)) / &
                      lai_f(i)**0.25_rkx/ustar(i,j)
            rcuto_f = max(rcuto_f, 100._rkx)
            rcuts_f = max(rcuts_f, 100._rkx)
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          end if
          !================================================================
          ! If snow occurs, Rg and Rcut are adjusted by snow cover
          ! fraction
          !================================================================
          fsnow = sd(i)/sdmax(kcov)
          fsnow = min(1.0_rkx, fsnow)   !snow cover fraction for leaves
!         print *, ' fsnow=  ', fsnow
          if ( fsnow > 0.0001_rkx .and. kcov /= 20 .or. &
                                        kcov /= 15 .or. &
                                        kcov /= 14 .or. &
                                        kcov /= 12 ) then
            rsnows = min(70.0_rkx*(275.15_rkx-ts(i)), 500._rkx)
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
            di = 0.001_rkx*ts(i)**1.75_rkx * &
              sqrt((29.0_rkx + mw(ig))/mw(ig)/29._rkx)
            di = di/1.0_rkx/(dair**0.3333_rkx + dgas**0.3333_rkx)**2
            vi = 145.8_rkx * 1.e-4_rkx * &
                 (ts(i) * 0.5_rkx + t2(i) *0.5_rkx)**1.5_rkx/ &
                 (ts(i) * 0.5_rkx + t2(i) *0.5_rkx + 110.4_rkx)
            !================================================================
            ! Calculate quasi-laminar resistance
            !================================================================
            rb(ig,i,j) = 5.0_rkx/ustar(i,j) * (vi/di)**.666667_rkx
!           print *, 'rb==', rb(ig,i,j)
            !================================================================
            ! Calculate stomatal resistance for each species from the ratio
            ! of  diffusity of water vapor to the gas species
            !================================================================
            dvh2o = 0.001_rkx*ts(i)**1.75_rkx * &
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
!           rc(ig,i,j) = (1.0_rkx - wst)/rstom + 1.0_rkx/(rg)+1.0_rkx/rcut
            rc(ig,i,j) = (1.0_rkx - wst)/rstom + 1.0_rkx/(rac+rg)+1.0_rkx/rcut
            rc(ig,i,j) = max(10._rkx,1.0_rkx/rc(ig,i,j))
          end do !igas
        end do !ilg
      end do !luc
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine stomtresis

end module mod_che_drydep
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
