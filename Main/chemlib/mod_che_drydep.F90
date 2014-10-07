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
  real(rk8) , parameter :: a1 = 145.8D0
  real(rk8) , parameter :: a2 = 1.5D0
  real(rk8) , parameter :: a3 = 110.4D0
!
! Molecular Free Path calculation parameters
!
  real(rk8) , parameter :: c1 = 6.54D-8
  real(rk8) , parameter :: c2 = 1.818D-5
  real(rk8) , parameter :: c3 = 1.013D5
  real(rk8) , parameter :: c4 = 293.15D0
!
! Cunningham slip correction factor parameters
!
  real(rk8) , parameter :: aa1 = 1.257D0
  real(rk8) , parameter :: aa2 = 0.4D0
  real(rk8) , parameter :: aa3 = 1.1D0
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
  real(rk8), parameter :: rainthr = 0.1D0
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
  real(rk8) lai(20,15)
  real(rk8) z01(20) , z02(20)
  integer(ik4) :: kk

  data (lai(1,kk), kk = 1, 15)/                    &
           0.1D0 , 0.1D0 , 0.1D0 , 0.5D0 , 1.0D0 , &
           2.0D0 , 3.0D0 , 3.5D0 , 4.0D0 , 0.1D0 , &
           0.1D0 , 0.1D0 , 0.1D0 , 0.1D0 , 4.0D0 /

  data (lai(2,kk), kk = 1, 15)/                    &
           1.0D0 , 1.0D0 , 1.0D0 , 1.0D0 , 1.0D0 , &
           1.0D0 , 1.0D0 , 1.0D0 , 1.0D0 , 1.0D0 , &
           1.0D0 , 1.0D0 , 1.0D0 , 1.0D0 , 1.0D0 /

  data (lai(3,kk), kk = 1, 15)/                    &
           5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , &
           5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , &
           5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 /

  data (lai(4,kk), kk = 1, 15)/                    &
           0.1D0 , 0.1D0 , 0.5D0 , 1.0D0 , 2.0D0 , &
           4.0D0 , 5.0D0 , 5.0D0 , 4.0D0 , 2.0D0 , &
           1.0D0 , 0.1D0 , 0.1D0 , 0.1D0 , 5.0D0 /

  data (lai(5,kk), kk = 1, 15)/                    &
           0.1D0 , 0.1D0 , 0.5D0 , 1.0D0 , 2.0D0 , &
           4.0D0 , 5.0D0 , 5.0D0 , 4.0D0 , 2.0D0 , &
           1.0D0 , 0.1D0 , 0.1D0 , 0.1D0 , 5.0D0 /

  data (lai(6,kk), kk = 1, 15)/                    &
           6.0D0 , 6.0D0 , 6.0D0 , 6.0D0 , 6.0D0 , &
           6.0D0 , 6.0D0 , 6.0D0 , 6.0D0 , 6.0D0 , &
           6.0D0 , 6.0D0 , 6.0D0 , 6.0D0 , 6.0D0 /

  data (lai(7,kk), kk = 1, 15)/                    &
           0.5D0 , 0.5D0 , 0.5D0 , 0.5D0 , 0.5D0 , &
           0.5D0 , 1.0D0 , 2.0D0 , 2.0D0 , 1.5D0 , &
           1.0D0 , 1.0D0 , 0.5D0 , 0.5D0 , 2.0D0  /

  data (lai(8,kk), kk = 1, 15)/                    &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 /

  data (lai(9,kk), kk = 1, 15)/                    &
           1.0D0 , 1.0D0 , 0.5D0 , 0.1D0 , 0.1D0 , &
           0.1D0 , 0.1D0 , 1.0D0 , 2.0D0 , 1.5D0 , &
           1.5D0 , 1.0D0 , 1.0D0 , 0.1D0 , 2.0D0 /

  data (lai(10,kk), kk = 1, 15)/                   &
           1.0D0 , 1.0D0 , 1.0D0 , 1.0D0 , 1.0D0 , &
           1.0D0 , 1.0D0 , 1.0D0 , 1.0D0 , 1.0D0 , &
           1.0D0 , 1.0D0 , 1.0D0 , 1.0D0 , 1.0D0 /

  data (lai(11,kk), kk = 1, 15)/                   &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 /

  data (lai(12,kk), kk = 1, 15)/                   &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 /

  data (lai(13,kk), kk = 1, 15)/                   &
           4.0D0 , 4.0D0 , 4.0D0 , 4.0D0 , 4.0D0 , &
           4.0D0 , 4.0D0 , 4.0D0 , 4.0D0 , 4.0D0 , &
           4.0D0 , 4.0D0 , 4.0D0 , 4.0D0 , 4.0D0 /

  data (lai(14,kk), kk = 1, 15)/                   &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 /

  data (lai(15,kk), kk = 1, 15)/                   &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 /

  data (lai(16,kk), kk = 1, 15)/                   &
           3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 , &
           3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 , &
           3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 /

  data (lai(17,kk), kk = 1, 15)/                   &
           3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 , &
           3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 , &
           3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 /

  data (lai(18,kk), kk = 1, 15)/                  &
           3.0D0 , 3.0D0 , 3.0D0 , 4.0D0 , 4.5D0 ,&
           5.0D0 , 5.0D0 , 5.0D0 , 4.0D0 , 3.0D0 ,&
           3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 , 5.0D0 /

  data (lai(19,kk), kk = 1, 15)/                  &
           3.0D0 , 3.0D0 , 3.0D0 , 4.0D0 , 4.5D0 ,&
           5.0D0 , 5.0D0 , 5.0D0 , 4.0D0 , 3.0D0 ,&
           3.0D0 , 3.0D0 , 3.0D0 , 3.0D0 , 5.0D0 /

  data (lai(20,kk), kk = 1, 15)/                  &
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 ,&
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 ,&
           0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 /

  data  z01/0.02D0, 0.04D0, 0.90D0, 0.40D0, 0.40D0, 2.00D0,&
            0.02D0, 0.04D0, 0.03D0, 0.05D0, 0.04D0, 0.01D0,&
            0.10D0, 0.00D0, 0.00D0, 0.20D0, 0.20D0, 0.60D0,&
            0.60D0, 0.00D0 /

  data  z02/0.10D0, 0.04D0, 0.90D0, 0.90D0, 1.00D0, 2.00D0,&
            0.10D0, 0.04D0, 0.03D0, 0.05D0, 0.04D0, 0.01D0,&
            0.10D0, 0.00D0, 0.00D0, 0.20D0, 0.20D0, 0.90D0,&
            0.90D0, 0.00D0 /
!
! Zhang stomatal resistance parameters
!
  real(rk8) :: tmin(20) , tmax(20)
  real(rk8) :: rsminz(20) , brs(20)
  real(rk8) :: topt(20) , bvpd(20)
  real(rk8) :: psi1(20) , psi2(20)
  real(rk8) :: rac1(20) , rac2(20)
  real(rk8) :: rgo(20) , rcutdO(20)
  real(rk8) :: rcutwO(20) , rcutdS(20)
  real(rk8) :: rgs(20) , sdmax(20)
  real(rk8) :: mw(31) , rm(31)
  real(rk8) :: alphaz(31) , betaz(31)

  data tmin /  5.0D0,    5.0D0,   -5.0D0,   -5.0D0,    0.0D0, &
               0.0D0,    5.0D0, -999.0D0,   -5.0D0,    5.0D0, &
            -999.0D0, -999.0D0,    0.0D0, -999.0D0, -999.0D0, &
               0.0D0,    0.0D0,   -3.0D0,    0.0D0, -999.0D0 /

  data tmax / 45.0D0,   40.0D0,   40.0D0,   40.0D0,   45.0D0, &
              45.0D0,   45.0D0, -999.0D0,   40.0D0,   45.0D0, &
            -999.0D0, -999.0D0,   45.0D0, -999.0D0, -999.0D0, &
              45.0D0,   45.0D0,   42.0D0,   45.0D0, -999.0D0 /

  data rsminz / 120.0D0, 150.0D0, 250.0D0, 250.0D0, 150.0D0, &
                150.0D0, 100.0D0,-999.0D0, 150.0D0, 150.0D0, &
               -999.0D0,-999.0D0, 150.0D0,-999.0D0,-999.0D0, &
                150.0D0, 250.0D0, 150.0D0, 150.0D0,-999.0D0 /

  data brs / 40.0D0,  50.0D0,  44.0D0,  44.0D0,  43.0D0, &
             40.0D0,  20.0D0,-999.0D0,  25.0D0,  40.0D0, &
           -999.0D0,-999.0D0,  40.0D0,-999.0D0,-999.0D0, &
             40.0D0,  44.0D0,  44.0D0,  43.0D0,-999.0D0 /

  data topt / 27.0D0,  30.0D0,  15.0D0,  15.0D0,  27.0D0, &
              30.0D0,  25.0D0,-999.0D0,  20.0D0,  25.0D0, &
            -999.0D0,-999.0D0,  20.0D0,-999.0D0,-999.0D0, &
              30.0D0,  25.0D0,  21.0D0,  25.0D0,-999.0D0 /

  data bvpd /  0.00D0,   0.00D0,   0.31D0,   0.31D0,   0.36D0, &
               0.27D0,   0.00D0,-999.00D0,   0.24D0,   0.09D0, &
            -999.00D0,-999.00D0,   0.27D0,-999.00D0,-999.00D0, &
               0.27D0,   0.27D0,   0.34D0,   0.31D0,-999.00D0 /

  data psi1 / -1.5D0,  -1.5D0,  -2.0D0,  -2.0D0,  -1.9D0, &
              -1.0D0,  -1.5D0,-999.0D0,   0.0D0,  -1.5D0, &
            -999.0D0,-999.0D0,  -1.5D0,-999.0D0,-999.0D0, &
              -2.0D0,  -2.0D0,  -2.0D0,  -2.0D0,-999.0D0 /

  data psi2 / -2.5D0,  -2.5D0,  -2.5D0,  -2.5D0,  -2.5D0, &
              -5.0D0,  -2.5D0,-999.0D0,  -1.5D0,  -2.5D0, &
            -999.0D0,-999.0D0,  -2.5D0,-999.0D0,-999.0D0, &
              -4.0D0,  -3.5D0,  -2.5D0,  -3.0D0,-999.0D0 /

  data rac1 / 10.0D0, 20.0D0,100.0D0, 60.0D0, 40.0D0, &
             250.0D0, 10.0D0,  0.0D0, 40.0D0, 20.0D0, &
               0.0D0,  0.0D0, 20.0D0,  0.0D0,  0.0D0, &
              60.0D0, 40.0D0,100.0D0,100.0D0,  0.0D0 /

  data rac2 / 40.0D0, 20.0D0,100.0D0,100.0D0, 40.0D0, &
             250.0D0, 40.0D0,  0.0D0, 40.0D0, 20.0D0, &
               0.0D0,  0.0D0, 20.0D0,  0.0D0,  0.0D0, &
              60.0D0, 40.0D0,100.0D0,100.0D0,  0.0D0 /

  data rcutdO / 4000.0D0, 4000.0D0, 4000.0D0, 4000.0D0, 6000.0D0, &
                6000.0D0, 4000.0D0, -999.0D0, 8000.0D0, 4000.0D0, &
                -999.0D0, -999.0D0, 5000.0D0, -999.0D0, -999.0D0, &
                6000.0D0, 5000.0D0, 4000.0D0, 4000.0D0, -999.0D0 /

  data rcutwO / 200.0D0, 200.0D0, 200.0D0, 200.0D0, 400.0D0, &
                400.0D0, 200.0D0,-999.0D0, 400.0D0, 200.0D0, &
               -999.0D0,-999.0D0, 300.0D0,-999.0D0,-999.0D0, &
                400.0D0, 300.0D0, 200.0D0, 200.0D0,-999.0D0 /

  data rgO / 200.0D0,  200.0D0,  200.0D0,  200.0D0,  200.0D0, &
             200.0D0,  200.0D0,  500.0D0,  500.0D0,  500.0D0, &
             500.0D0, 2000.0D0,  500.0D0, 2000.0D0, 2000.0D0, &
             200.0D0,  200.0D0,  200.0D0,  200.0D0, 2000.0D0 /

  data rcutds / 1500.0D0, 1000.0D0, 2000.0D0, 2000.0D0, 2500.0D0, &
                2500.0D0, 1000.0D0, -999.0D0, 2000.0D0, 2000.0D0, &
                -999.0D0, -999.0D0, 1500.0D0, -999.0D0, -999.0D0, &
                2000.0D0, 2000.0D0, 2500.0D0, 2500.0D0, -999.0D0 /

  data rgs / 200.0D0, 200.0D0, 200.0D0, 200.0D0, 200.0D0, &
             100.0D0, 200.0D0, 700.0D0, 300.0D0,  50.0D0, &
             700.0D0,  70.0D0,  50.0D0,  20.0D0,  20.0D0, &
             200.0D0, 200.0D0, 200.0D0, 200.0D0,  20.0D0 /

  data sdmax /  10.0D0,  5.0D0, 200.0D0,   1.1D0, 200.0D0, &
               400.0D0, 20.0D0,   2.0D0,   2.0D0,  10.0D0, &
                 2.0D0,  1.0D0,  10.0D0,-999.0D0,-999.0D0, &
                50.0D0, 50.0D0, 200.0D0, 200.0D0,-999.0D0 /

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

  data rm / 0.0D0 ,   0.0D0 ,   0.0D0 ,   0.0D0 ,    0.0D0 , &
            0.0D0 ,   0.0D0 ,   0.0D0 ,   0.0D0 ,    0.0D0 , &
            0.0D0 ,   0.0D0 ,   0.0D0 ,   0.0D0 ,  100.0D0 , &
          100.0D0 , 100.0D0 , 100.0D0 , 100.0D0 ,    0.0D0 , &
          100.0D0 ,   0.0D0   , 0.0D0 ,   0.0D0 ,    0.0D0 , &
            0.0D0 ,   0.0D0   , 0.0D0 ,   0.0D0 ,  100.0D0 , &
          100.0D0 /

  data alphaz /  1.00D0 , 1.00D0 , 0.00D0 , 0.00D0 , 1.00D0 , &
                10.00D0 , 2.00D0 , 5.00D0 , 1.00D0 , 0.00D0 , &
                 0.00D0 , 0.00D0 , 0.00D0 , 0.80D0 , 0.00D0 , &
                 0.00D0 , 0.00D0 , 0.00D0 , 0.00D0 , 0.00D0 , &
                 0.00D0 , 0.01D0 , 0.60D0 , 0.60D0 , 0.40D0 , &
                 0.01D0 , 2.00D0 , 1.50D0 , 0.10D0 , 0.00D0 , &
                 0.00D0 /

  data betaz /  0.00D0 , 1.00D0 , 0.80D0 , 1.00D0 , 1.00D0 , &
               10.00D0 , 2.00D0 , 5.00D0 , 0.00D0 , 0.60D0 , &
                0.60D0 , 0.80D0 , 0.30D0 , 0.20D0 , 0.05D0 , &
                0.05D0 , 0.05D0 , 0.05D0 , 0.05D0 , 0.05D0 , &
                0.05D0 , 0.00D0 , 0.10D0 , 0.00D0 , 0.00D0 , &
                0.00D0 , 0.00D0 , 0.00D0 , 0.80D0 , 0.50D0 , &
                0.50D0 /


  data mw / 64.0D0 ,  98.0D0 ,  46.0D0 ,  48.0D0 ,  34.0D0 , &
            63.0D0 ,  47.0D0 ,  79.0D0 ,  17.0D0 , 121.0D0 , &
           135.0D0 , 183.0D0 , 147.0D0 ,  30.0D0 ,  44.0D0 , &
            58.0D0 ,  72.0D0 , 128.0D0 , 106.0D0 ,  70.0D0 , &
            70.0D0 ,  72.0D0 ,  32.0D0 ,  46.0D0 ,  60.0D0 , &
           104.0D0 ,  46.0D0 ,  60.0D0 ,  48.0D0 ,  77.0D0 , &
           147.0D0 /

  contains

    subroutine drydep_aero(j,mbin,indsp,rhop,ivegcov,throw,roarow, &
                           shj,pressg,temp2,sutemp,srad,rh10,      &
                           wind10,zeff,beffdiam,pdepv,ddepv)
!
      implicit none
!
      integer(ik4) , intent(in) :: j , mbin
      integer(ik4) , intent(in) , dimension(mbin) :: indsp
      integer(ik4) , intent(in) , dimension(ici1:ici2) :: ivegcov
      real(rk8) , dimension(ici1:ici2) , intent(in) :: pressg , rh10 , &
                       srad , sutemp , temp2 , wind10 , zeff
      real(rk8) , dimension(ici1:ici2,kz) , intent(in) :: roarow , throw
      real(rk8) , dimension(kz) , intent(in) :: shj
      real(rk8) , dimension(mbin) , intent(in) :: beffdiam
      real(rk8) , intent(in) :: rhop

! output table to be passed out. Care dimension is ntr

      real(rk8) , intent(out) , dimension(ici1:ici2,kz,ntr) :: pdepv
      real(rk8) , intent(out) , dimension(ici1:ici2,ntr) :: ddepv
!
      real(rk8) :: amfp , amob , eb , eim , ein , frx1
      real(rk8) :: pre , prii , priiv , r1 , st
      real(rk8) , dimension(ici1:ici2,kz) :: amu
      real(rk8) , dimension(ici1:ici2) :: anu , schm
      real(rk8) , dimension(ici1:ici2,kz,mbin) :: cfac , pdepvsub , pdiff , &
                  rhsize , taurel
      real(rk8) , dimension(ici1:ici2,luc) :: ra , ustar
      real(rk8) , dimension(ici1:ici2,luc,mbin) :: rs
      real(rk8), dimension(ici1:ici2,kz) :: wk, settend
      real(rk8) , dimension(mbin) :: avesize
      integer(ik4) :: i , k , kcov , l , n , ib
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'drydep_aero'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif
      ! here avesize is a RADIUS of the particle bin in m :
      ! calculated here from bin effective diameter in micrometer
      do n = 1 , mbin
        avesize(n) = beffdiam(n)* 1.D-06 * d_half
      end do
!
!======================================================================
!
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
            amu(i,l) = a1*1.D-8*throw(i,l)**a2/(throw(i,l)+a3)
            ! mid layer pressure in [pascal].
            pre = pressg(i)*shj(l)
            !
            ! ********************************************************
            ! * mean molecular free path.                         ****
            ! *     k.v. beard [1976], j atm. sci., 33            ****
            ! ********************************************************
            !
            amfp = c1*(amu(i,l)/c2)*(c3/pre)*sqrt(throw(i,l)/c4)
            prii = 2.0D0/9.0D0*egrav/amu(i,l)
            priiv = prii*(rhop-roarow(i,l))
            !
            ! ********************************************************
            ! * cunningham slip correction factor and             ****
            ! * relaxation time = vg/grav.                        ****
            ! ********************************************************
            !
            cfac(i,l,n) = 1.0D0 + amfp/avesize(n) * &
                         (aa1+aa2*exp(-aa3*avesize(n)/amfp))
            taurel(i,l,n) = dmax1(priiv*avesize(n)**2*cfac(i,l,n) * &
                            regrav,0.D0)
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
            frx1 = 1.0D0
            rhsize(i,l,n) = avesize(n)*frx1 ! still a radius
            anu(i) = amu(i,l)/roarow(i,l)
            amob = 6.0D0*mathpi*amu(i,l)*rhsize(i,l,n)/cfac(i,l,n)
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
                eb = schm(i)**(-0.666667D0)
!               eim=(st/(st+aest(k)))**2
                eim = (st/(st+aest(kcov)))**2
                eim = dmin1(eim,0.6D0)
                ein = 0.0D0
!               if (arye(k) > 0.0001D0) then
!                 ein = (1000.0D0*2.0D0*avesize(n)/arye(k))**1.5D0
!               end if
                if ( arye(kcov) > 0.0001D0 ) then
                  ein = (1000.0D0*2.0D0*avesize(n)/arye(kcov))**1.5D0
                end if
                ein = dmin1(ein,0.5D0)
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

                r1 = dmax1(0.5D0,exp(-(min(sqrt(st),25.0D0))))
                if ( kcov >= 11 .and. r1 < 0.5D0 ) r1 = 0.5D0
                if ( r1 < 0.4D0 ) r1 = 0.4D0
                ! ***************************************************
                ! * calculation of rs: the surface resistance   *****
                ! * which depends on the collection efficiency  *****
                ! * of the surface and is determined by the     *****
                ! * various deposition processes                *****
                ! ***************************************************
!               rs= 1.0/ustar(i,k)/(eb+eim+ein)/r1
                rs(i,k,n) = 1.0D0/3.0D0/ustar(i,k)/(eb+eim+ein)/r1
              end do
            end do
          end if
        end do
      end do

!======================================================================


      ! average settling and deposition velocities on bin
      ! care we use pdepv and ddpv table that are dimensionned to ntr
      ! and not mbin !
      do ib = 1 , mbin
        ! there isw no sub-bin anymore / we consider directly effective radius
        pdepv(:,:,indsp(ib)) = 0.0D0
        ddepv(:,indsp(ib))   = 0.0D0
        do i = ici1 , ici2
          pdepv(i,:,indsp(ib)) = pdepvsub(i,:,ib)
          ! agregate the dry deposition velocity, remember one cover per grid
          ! cell for now
          ! the dry deposition deposition velocity must accound also for the
          ! settling vrlocity at kz
          ! simple form now add the vs
          ddepv(i,indsp(ib)) = 1.0D0/(ra(i,1)+rs(i,1,ib)) + pdepvsub(i,kz,ib)
        end do
      end do
      !
      ! Finally update the emission and settling tendencies for
      ! dust and sea salt
      !
      do ib = 1 , mbin
        ! deposition, remember chiten must be normalised by psb and consistent with chib
        do k = 2 , kz
          do i = ici1 , ici2
            wk(i,k) =  twt(k,1)*chib(j,i,k,indsp(ib))+ &
                       twt(k,2)*chib(j,i,k-1,indsp(ib))
          end do
        end do
        do i = ici1 , ici2
          do k = 2 , kz - 1
            ! do not apply to the first level
            settend(i,k) = (wk(i,k+1)*pdepv(i,k+1,indsp(ib)) - &
                            wk(i,k)*pdepv(i,k,indsp(ib))) / cdzq(j,i,k)

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

            settend(i,kz) =  (chib(j,i,kz,indsp(ib))  * ddepv(i,indsp(ib))-  &
                              wk(i,kz)*pdepv(i,kz,indsp(ib))) / cdzq(j,i,kz)

            chiten(j,i,kz,indsp(ib)) = chiten(j,i,kz,indsp(ib)) - settend(i,kz)

            !diagnostic for settling and drydeposition removal
            if ( ichdiag == 1 ) then
              cseddpdiag(j,i,kz,indsp(ib)) = cseddpdiag(j,i,kz,indsp(ib)) - &
                                             settend(i,kz) * cfdout
            end if
            ! diagnostic for dry deposition flux (in kg .m2.s-1) accumulated
            remdrd(j,i,indsp(ib)) = remdrd(j,i,indsp(ib)) + &
                 chib3d(j,i,kz,indsp(ib)) *crhob3d(j,i,kz)* &
                 ddepv(i,indsp(ib)) * cfdout
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
!
!
!******************************************************************************
!******************************************************************************
!
    subroutine drydep_gas(j,lmonth,lday,ivegcov,rh10,srad,tsurf, &
                          prec,temp10,wind10,zeff)

      use mod_che_indices
      implicit none
      integer(ik4) , intent(in) :: j
      integer, intent(in) :: lmonth , lday
      integer(ik4) , intent(in) , dimension(ici1:ici2) :: ivegcov
      real(rk8) , intent(in) , dimension(ici1:ici2) :: rh10 , srad , tsurf , &
                                            prec, temp10 , wind10 , zeff
      real(rk8),  dimension(ici1:ici2,ntr) :: drydepvg

      integer(ik4) :: n , i , im , kcov
      real(rk8) , dimension(ici1:ici2,luc) :: ustar, resa
      real(rk8) , dimension(ngasd,ici1:ici2,luc) :: resb, resc
      real(rk8) , dimension(ngasd,ici1:ici2,luc) :: vdg
      real(rk8) , dimension(ici1:ici2) :: icz , ddrem
      real(rk8) , dimension(ici1:ici2) :: lai_f , laimin , laimax , snow
      real(rk8) :: kd
#ifdef CLM
      integer(ik4) :: jj, ii
#endif
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
                     lai(kcov,im))/30.D0 * dble(15 + lday)
        else
          lai_f(i) = lai(kcov,lmonth) + (lai(kcov,lmonth+1) - &
                     lai(kcov,lmonth))/30.D0 * dble(lday - 15)
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
!         drydepvg(:,inh3)  =  vdg(9,:,1)!*0.5
         drydepvg(:,ipan)  =  vdg(10,:,1)!*0.5
         drydepvg(:,ihcho) =  vdg(14,:,1)!*0.5
         drydepvg(:,iald2) =  vdg(15,:,1)!*0.5
         drydepvg(:,ich3oh)  =  vdg(23,:,1)!*0.5
       end if

       ! Finally : gas phase dry dep tendency calculation
#ifdef CLM
       jj = global_cross_jstart+j-1
#endif
       if ( ichdrdepo == 1 ) then
         do i = ici1 , ici2

           ! If using CLM then use the dry deposition velocities coming directly
           ! from internal CLM calculations
#ifdef CLM
           if ( ivegcov(i) == 0 ) then
              ii = global_cross_istart+i-1
              drydepvg(i,:)  =  cdep_vels(jj,ii,:)
           end if
#endif
           do n = 1 , ntr
             kd =  drydepvg(i,n) / cdzq(j,i,kz) !Kd removal rate in s-1
             if ( kd*dt < 25.0D0 ) then
               ! dry dep removal tendency (+)
               ddrem(i) = chib(j,i,kz,n) * (d_one-dexp(-kd*dt))/dt
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
             drydepv(j,i,n) =  d_zero
             ddv_out(j,i,n) =  ddv_out(j,i,n) + drydepvg(i,n)
           end do
         end do
       else if ( ichdrdepo == 2 ) then
         do i = ici1 , ici2

           ! If using CLM then use the dry deposition velocities coming directly
           ! from internal CLM calculations
#ifdef CLM
           if ( ivegcov(i) == 0 ) then
              ii = global_cross_istart+i-1
              drydepvg(i,:)  =  cdep_vels(jj,ii,:)
           end if
#endif
           do n = 1 , ntr
             chifxuw(j,i,n) = chifxuw(j,i,n) - chib(j,i,kz,n) / &
                               cpsb(j,i) * drydepvg(i,n)
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

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    subroutine aerodyresis(zeff,wind10,temp2,sutemp,rh10,srad,ivegcov,ustar,ra)
      implicit none
      integer(ik4) , dimension(ici1:ici2) , intent(in) :: ivegcov
      real(rk8) , dimension(ici1:ici2) , intent(in) :: temp2 , wind10 , rh10
      real(rk8) , dimension(ici1:ici2) , intent(in) :: sutemp , srad , zeff
      real(rk8) , dimension(ici1:ici2,luc) , intent(out) :: ustar , ra
!
      integer(ik4) :: i , j
      real(rk8) :: vp , tsv
      real(rk8) :: z , zl , ww
      real(rk8) :: ptemp2 , es , qs
      real(rk8) :: wvpm , vptemp , tsw , mol
      real(rk8) :: z0water , dthv , cun , zdl
      real(rk8) :: psiu , psit , x , y
      real(rk8) :: thstar , rib , dtemp , tbar
      real(rk8) :: ustarsq , utstar , kui
      real(rk8) :: ratioz , logratio , asq
      real(rk8) :: aa , cm , ch , fm , fh
      real(rk8) , dimension(ici1:ici2) :: zz0
      real(rk8) , parameter :: z10 = 10.0D0
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'aerodyresis'
      integer(ik4) , save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif

      !======================================================================
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
          ww = dmax1(wind10(i),1.0D0)
          zz0(i) = zeff(i)
          ! ***************************************************************
          ! * potential temperature at z2  (deg. k)
          ! ***************************************************************
          ptemp2 = temp2(i) + z10*0.0098D0
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
            es = 6.108D0*exp(17.27D0*(temp2(i)-tzero)/(temp2(i)-35.86D0))
            vp = rh10(i)*es
            wvpm = ep2*vp/(stdpmb-vp)
            vptemp = ptemp2*(1.0D0+0.61D0*wvpm)
            ! **************************************************************
            ! * assume rh10 at water surface is 100%
            ! *   vp = es(tsw-tzero) !sat. vap press at surface
            ! *   saturated vapour pressure at surface
            ! *   saturated mixing ratio at surface
            ! *   tsv - virtual potential temperature at surface (deg. k)
            ! **************************************************************
            tsw = sutemp(i)
            vp = 6.108D0*exp(17.27D0*(tsw-tzero)/(tsw-35.86D0))
            qs = ep2*vp/(stdpmb-vp)
            tsv = tsw*(1.0D0+0.61D0*qs)
            z0water = 1.0D-4
            ! **************************************************************
            ! * scalet  :  not required if  z2 = 10m
            ! **************************************************************
            dthv = (vptemp-tsv)
            ! **************************************************************
            ! * calculate drag coefficient cun with neutral condition
            ! * assumption  garratt (1977)
            ! **************************************************************
            cun = 7.5D-4 + 6.7D-5*ww
            mol = 9999.0D0
            if ( abs(dthv) > 1.0D-6 ) then
              mol = vptemp*cun**1.5D0*ww**2/(5.096D-3*dthv)
            end if
            if ( mol > 0.0D0  .and. mol < 5.0D0 ) mol =  5.0D0
            if ( mol > -5.0D0 .and. mol < 0.0D0 ) mol = -5.0D0
            zdl = z10/mol
            if ( zdl < 0.0D0 ) then
              ! **************************************************************
              ! * wind speed
              ! **************************************************************
              x = (1.0D0-15.0D0*zdl)**0.25D0
              psiu = 2.0D0*dlog(0.5D0*(1.0D0+x))+dlog(0.5D0*(1.0D0+x*x)) - &
                     2.0D0*atan(x) + 0.5D0*mathpi
              ! **************************************************************
              ! * pot temp
              ! **************************************************************
              y = sqrt(1.0D0-9.0D0*zdl)
              psit = 2.0D0*0.74D0*dlog((1.0D0+y)/2.0D0)
            else
              psiu = -4.7D0*zdl
              psit = psiu
            end if
            z0water = 0.000002D0*ww**2.5D0
            ustar(i,j) = vonkar*ww/(dlog(z10/z0water)-psiu)
            thstar = vonkar*(ptemp2-sutemp(i)) / &
                     (0.74D0*dlog(z10/z0water)-psit)
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
            if ( srad(i) > 0.0D0 .and. rib > 0.0D0 ) rib = 1.D-15
            dtemp = ptemp2 - sutemp(i)
            if ( dabs(dtemp) < 1.D-10 ) dtemp = dsign(1.D-10,dtemp)
            tbar = 0.5D0*(ptemp2+sutemp(i))
            ratioz = z10/zz0(i)
            logratio = dlog(ratioz)
            asq = 0.16D0/(logratio**2)
            if ( rib <= 0.0D0 ) then
              aa = asq*9.4D0*dsqrt(ratioz)
              cm = 7.4D0*aa
              ch = 5.3D0*aa
              fm = 1.0D0 - (9.4D0*rib/(1.0D0+cm*dsqrt(dabs(rib))))
              fh = 1.0D0 - (9.4D0*rib/(1.0D0+ch*dsqrt(dabs(rib))))
            else
              fm = 1.0D0/((1.0D0+4.7D0*rib)**2)
              fh = fm
            end if
            ustarsq = asq*ww**2*fm
            utstar = asq*ww*dtemp*fh/0.74D0
            ustar(i,j) = dsqrt(ustarsq)
            thstar = utstar/ustar(i,j)
            mol = tbar*ustarsq/(vonkar*egrav*thstar)
          end if

          kui = 1.0D0/(vonkar*ustar(i,j))

          ! **************************************************************
          ! * compute the values of  ra                            *******
          ! **************************************************************
          z = z10
          zl = z/mol
          if ( zl >= 0.0D0 ) then
            ra(i,j) = kui*(0.74D0*dlog(z/zz0(i))+4.7D0*zl)
          else
            ra(i,j) = kui*0.74D0*(dlog(z/zz0(i))- &
                      2.0D0*dlog((1.0D0+sqrt(1.0D0-9.0D0*zl))*0.5D0))
          end if
          ra(i,j) = dmax1(ra(i,j),0.99D0)
          ra(i,j) = dmin1(ra(i,j),999.9D0)
        end do
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine aerodyresis

!***********************************************************************

    subroutine stomtresis(lai_f,laimin,laimax,ivegcov,igas, &
                          ustar,prec,sd,srad,ts,t2,rh,coszen,rc,rb)

      implicit none

      integer(ik4) , intent(in) :: igas
      integer(ik4) , intent(in) , dimension(ici1:ici2) :: ivegcov
      real(rk8) , dimension(ici1:ici2) , intent(in) :: coszen, srad , ts , rh , &
                                               prec , sd , t2
      real(rk8) , dimension(ici1:ici2) , intent(in) :: lai_f , laimin , laimax
      real(rk8) , intent(in) , dimension(ici1:ici2,luc) :: ustar
      real(rk8) , intent(out) , dimension(igas,ici1:ici2,luc) :: rb , rc
!
      integer(ik4) :: i , j , kcov , ig
      real(rk8) :: rst , wst , rac , rgs_f
      real(rk8) :: rdu , rdv , rgo_f
      real(rk8) :: rcuto_f , rcuts_f
      real(rk8) :: ww1 , ww2 , ww3
      real(rk8) :: rdm , rdn , rv , rn
      real(rk8) :: ratio , sv , fv , fvv
      real(rk8) :: pardir , pardif
      real(rk8) :: tmaxk , tmink
      real(rk8) :: pshad , psun , rshad , rsun
      real(rk8) :: gshad , gsun , fsun , fshad
      real(rk8) :: gspar , temps !C
      real(rk8) :: bt , gt , gw , ryx
      real(rk8) :: es , d0 , gd , psi
      real(rk8) :: coedew , dq , usmin
      real(rk8) :: fsnow , rsnows
      real(rk8) :: dgas , di , vi
      real(rk8) :: dvh2o , rstom
      real(rk8) :: rcut , rg , xp
      logical :: is_dew , is_rain
      real(rk8) , parameter :: dair = 0.369D0 * 29.0D0 + 6.29D0
      real(rk8) , parameter :: dh2o = 0.369D0 * 18.0D0 + 6.29D0
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
          rst = -999.0D0
          if (srad(i)   >= 0.1D0  .and. &
              ts(i)     <  tmaxk  .and. &
              ts(i)     >  tmink  .and. &
              lai_f(i)  > 0.001D0 .and. &
              coszen(i) > 0.001D0 ) then
            !================================================================
            ! Calculate direct and diffuse PAR from solar radiation and
            ! solar zenith angle
            !================================================================
            rdu   = 600.0D0 * dexp(-0.185D0/coszen(i))*coszen(i)
            rdv   = 0.4D0 * (600.0D0 - rdu ) * coszen(i)
            ww1   = -dlog(coszen(i))/2.302585D0
!           print *, ' ww1 = ', ww1
            ww2   = -1.195D0 + 0.4459D0 * ww1 - 0.0345D0 * ww1**2
            ww3   = 1320.0D0*10.0D0**ww2
!           print *, 'ww= ', ww
            rdm   = (720.0D0*dexp(-0.06D0/coszen(i))-ww3)*coszen(i)
!           print *, 'ww3= ', ww3, rdm
            rdn   = 0.6D0 * (720.0D0 - rdm - ww3) * coszen(i)
            rv    = dmax1(0.1D0,  rdu + rdv)
            rn    = dmax1(0.01D0, rdm + rdn)
            ratio = dmin1(0.9D0,srad(i)/( rv + rn))
!           print *, 'ratio= ', ratio, rdn, rv, rn
            sv    = ratio * rv                            ! Total PAR
            fv    = dmin1(0.99D0, (0.901D0 - ratio)/0.7D0)
!           print *, 'sv  fv  = ', sv, fv
!           print *, 'rv  xxxxx  = ', rv, (1.0 - fv**0.6667)
            fvv   = dmax1(0.01D0,rdu/rv*(1.0D0 - fv**0.6667D0))
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
            if ( lai_f(i) > 2.5D0 .and. srad(i) > 200.0D0 ) then
              pshad = pardif * dexp(-0.5D0 * lai_f(i)**0.8D0) + &
                      0.07D0 * pardir * (1.1D0-0.1D0*lai_f(i))* &
                      dexp(-coszen(i))
              psun = pardir**0.8D0*0.5D0/coszen(i) + pshad
            else
              pshad = pardif * dexp(-0.5D0 * lai_f(i)**0.7D0) + &
                      0.07D0 * pardir *(1.1D0-0.1D0*lai_f(i)) * &
                      dexp(-coszen(i))
              psun = pardir * 0.5D0/coszen(i) + pshad
            end if
!           print *, 'pshad   psun   ', pshad , psun
            rshad = rsminz(kcov) + brs(kcov) * rsminz(kcov)/pshad
            rsun  = rsminz(kcov) + brs(kcov) * rsminz(kcov)/psun
            gshad = 1.0D0/rshad
            gsun  = 1.0D0/rsun
!           print *, 'rshad  ----< ', rshad, rsun, ' >---------rsun'
!           print *, 'gshad  ----< ', gshad, gsun, ' >-------- gsun'
            !================================================================
            ! Fsun, Fshade are the total sunlit and shaded leaf area
            ! index
            !================================================================
            xp = 0.5D0*lai_f(i)/coszen(i)
            if ( xp < 25.0D0 ) then
              fsun  = 2.0D0*coszen(i)*(1.0D0-dexp(-xp))
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
            es = 6.108D0*dexp(17.27D0*(ts(i)-tzero)/(ts(i)-35.86D0))
            d0 = es*(d_one-rh(i))/10.0D0 ! kPa
            gd = 1.0D0 - bvpd(kcov) * d0
!           print *, 'gd===',gd
            !================================================================
            ! function for water stress
            !================================================================
            psi = (-0.72D0 - 0.0013D0 * srad(i))
!           psi_s = (-0.395D0-0.043D0*(ts-tzero))*102.0D0
            gw = (psi - psi2(kcov))/(psi1(kcov) - psi2(kcov))
!           print *, 'gw==',gw
!           TEST
!           gw = 1
            if ( gw > 1.0D0 ) gw = 1.0D0
            if ( gw < 0.1D0 ) gw = 0.1D0
            if ( gd > 1.0D0 ) gd = 1.0D0
            if ( gd < 0.1D0 ) gd = 0.1D0
            !================================================================
            ! Stomatal resistance for water vapor
            !================================================================
            rst = 1.0D0 / (gspar * gt * gd * gw)
!           print *, 'rst===',rst
          end if
          coedew = 0.1D0  ! for clear cloud
          es = 6.108D0*dexp(17.27D0*(ts(i)-tzero)/(ts(i)-35.86D0))
          dq = 0.622D0/1000.0D0*es*(1.0D0-rh(i))*1000.0D0 ! unit g/kg
          dq = dmax1(0.0001D0,dq)
          usmin = 1.5D0/dq*coedew
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
          wst = 0.0D0
          if ( (is_dew .or. is_rain) .and. srad(i) > 200.0D0 ) then
            wst = (srad(i) - 200.0D0)/800.0D0
            wst = dmin1(wst, 0.5D0)
          end if
          !================================================================
          ! In-canopy aerodynamic resistance
          !================================================================
          rac = rac1(kcov)+(lai_f(i)-laimin(i))/ &
                (laimax(i)-laimin(i)+1.D-10)*(rac2(kcov)-rac1(kcov))
!         print *, 'rac1 = ', rac
          rac = rac*lai_f(i)**0.25D0/ustar(i,j)/ustar(i,j)
!         print *, 'rac2 = ', rac
          !================================================================
          ! Ground resistance for O3
          !================================================================
          if (ts(i) < 272.15D0 .and. kcov /= 14 ) then
            rgo_f = dmin1( rgo(kcov)*2.0D0, rgo(kcov) *     &
                           dexp(0.2D0*(272.15D0-ts(i))))
!           print *, 'rgo_f1 =',rgo_f, ts(i)
          else
            rgo_f = rgo(kcov)
          end if
          !================================================================
          ! Ground resistance for SO2
          !================================================================
          if ( kcov == 12 ) then
            rgs_f = dmin1(rgs(kcov)*(275.15D0 - ts(i)), 500.D0)
            rgs_f = dmax1(rgs(kcov), 100.D0)
!           print *, 'rgs_f ==== ', rgs_f
          else if ( is_rain .and. kcov /= 14 ) then
            rgs_f = 50.0D0
!           print *, 'rgs_f ==== ', rgs_f
          else if ( is_dew .and. kcov /= 14 ) then
            rgs_f = 100.0D0
!           print *, 'rgs_f ==== ', rgs_f
          else if ( ts(i) < 272.156D0 .and. kcov /= 14 ) then
            rgs_f = dmin1(rgs(kcov)*2.0D0, rgs(kcov) *     &
                          dexp(0.2D0*(272.156D0 - ts(i))))
!           print *, 'rgs_f ==== ', rgs_f
          else
            rgs_f = rgs(kcov)
!           print *, 'rgs_f ==== ', rgs_f
          end if
          !================================================================
          ! Cuticle resistance for O3 AND SO2
          !================================================================
          if ( rcutdo(kcov) <= -1.0D0 ) then
            rcuto_f = 1.D25
            rcuts_f = 1.D25
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          else if ( is_rain ) then
            rcuto_f = rcutwo(kcov)/sqrt(lai_f(i))/ustar(i,j)
            rcuts_f = 50.0D0/sqrt(lai_f(i))/ustar(i,j)
            rcuts_f = dmax1(rcuts_f, 20.D0)
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          else if ( is_dew ) then
            rcuto_f = rcutwo(kcov)/sqrt(lai_f(i))/ustar(i,j)
            rcuts_f = 100.0D0/sqrt(lai_f(i))/ustar(i,j)
            rcuts_f = dmax1(rcuts_f, 20.D0)
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          else if (ts(i) < 272.156D0 ) then
            ryx = dexp(0.2D0 * (272.156D0 - ts(i) ))
            rcuto_f = rcutdo(kcov)/dexp(3.0D0 * rh(i))/     &
                      lai_f(i)**0.25D0/ustar(i,j)
            rcuts_f = rcutds(kcov)/dexp(3.0D0 * rh(i))/     &
                      lai_f(i)**0.25D0/ustar(i,j)
            rcuto_f = dmin1(rcuto_f * 2.0D0, rcuto_f * ryx )
            rcuts_f = dmin1(rcuts_f * 2.0D0, rcuts_f * ryx )
            rcuto_f = dmax1(rcuto_f,100.D0)
            rcuts_f = dmax1(rcuts_f,100.D0)
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          else
            rcuto_f = rcutdo(kcov)/exp(3.0D0*rh(i)) / &
                      lai_f(i)**0.25D0/ustar(i,j)
            rcuts_f = rcutds(kcov)/exp(3.0D0*rh(i)) / &
                      lai_f(i)**0.25D0/ustar(i,j)
            rcuto_f = dmax1(rcuto_f, 100.D0)
            rcuts_f = dmax1(rcuts_f, 100.D0)
!           print *, 'RCUT === ', rcuto_f,rcuts_f
          end if
          !================================================================
          ! If snow occurs, Rg and Rcut are adjusted by snow cover
          ! fraction
          !================================================================
          fsnow = sd(i)/sdmax(kcov)
          fsnow = dmin1(1.0D0, fsnow)   !snow cover fraction for leaves
!         print *, ' fsnow=  ', fsnow
          if ( fsnow > 0.0001D0 .and. kcov /= 20 .or. &
                                      kcov /= 15 .or. &
                                      kcov /= 14 .or. &
                                      kcov /= 12 ) then
            rsnows = dmin1(70.0D0*(275.15D0-ts(i)), 500.D0)
            rsnows = dmax1(rsnows, 100.D0)
            rcuts_f = 1.0D0/((1.0D0 - fsnow)/rcuts_f + fsnow/rsnows)
            rcuto_f = 1.0D0/((1.0D0 - fsnow)/rcuto_f + fsnow/2000.0D0)
            fsnow = dmin1(1.0D0, fsnow*2.0D0)
            ! snow cover fraction for ground
            rgs_f = 1.0D0/((1.0D0 - fsnow)/rgs_f + fsnow/rsnows)
            rgo_f = 1.0D0/((1.0D0 - fsnow)/rgo_f + fsnow/2000.0D0)
!           print *, 'rsnows= ', rsnows, ' fsnow=  ', fsnow
          end if
          !================================================================
          ! Calculate diffusivity for each gas species
          !================================================================
          do ig = 1 , igas
            dgas = 0.369D0 * mw(ig) + 6.29D0
            di = 0.001D0*ts(i)**1.75D0*sqrt((29.0D0 + mw(ig))/mw(ig)/29.D0)
            di = di/1.0D0/(dair**0.3333D0 + dgas**0.3333D0)**2
            vi = 145.8D0 * 1.D-4 * (ts(i) * 0.5D0 + t2(i) *0.5D0)**1.5D0/ &
                 (ts(i) * 0.5D0 + t2(i) *0.5D0 + 110.4D0)
            !================================================================
            ! Calculate quasi-laminar resistance
            !================================================================
            rb(ig,i,j) = 5.0D0/ustar(i,j) * (vi/di)**.666667D0
!           print *, 'rb==', rb(ig,i,j)
            !================================================================
            ! Calculate stomatal resistance for each species from the ratio
            ! of  diffusity of water vapor to the gas species
            !================================================================
            dvh2o = 0.001D0*ts(i)**1.75D0*sqrt((29.0D0+18.0D0)/29.0D0/18.0D0)
            dvh2o = dvh2o/(dair**0.3333D0 + dh2o**0.3333D0)**2
            rstom = rst * dVh2o/di + rm(ig)
            ! (rst <999) for bare surfaces)
            !================================================================
            ! Scale cuticle and ground resistances for each species
            !================================================================
            rcut = 1.0D0/(alphaz(ig)/rcuts_f+betaz(ig)/rcuto_f)
            rg   = 1.0D0/(alphaz(ig)/rgs_f+betaz(ig)/rgo_f)
            !================================================================
            ! Calculate total surface resistance
            !================================================================
            ! account for zero stomatal resistance (rst and rstom are zero
            ! for bare surfaces)
            ! set wst to 1 also in that case (total stomatal blocking).
            if ( rst == -999.0 ) wst = 1.0D0
!           rc(ig,i,j) = (1.0D0 - wst)/rstom + 1.0D0/(rg)+1.0D0/rcut
            rc(ig,i,j) = (1.0D0 - wst)/rstom + 1.0D0/(rac+rg)+1.0D0/rcut
            rc(ig,i,j) = dmax1(10.D0,1.0D0/rc(ig,i,j))
          end do !igas
        end do !ilg
      end do !luc
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine stomtresis

end module mod_che_drydep
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
