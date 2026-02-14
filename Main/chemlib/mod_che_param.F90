!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_che_param

  use mod_intkinds
  use mod_realkinds

  implicit none

  public

  integer(ik4) :: ichaer
  !
  !FAB : redefine soil prop. similar to bats for chemistry externalisation.
  ! think about an interface!
  !
  !     ******      xmopor is fraction of soil that is voids
  real(rkx), dimension(12), parameter :: cxmopor = &
     [ 0.33_rkx, 0.36_rkx, 0.39_rkx, 0.42_rkx, 0.45_rkx, &
       0.48_rkx, 0.51_rkx, 0.54_rkx, 0.57_rkx, 0.60_rkx, &
       0.63_rkx, 0.66_rkx]

  integer(ik4), dimension(22), parameter :: ciexsol = &
     [ 6, 6, 6, 6, 7, 8, 6, 1, 6, 6, 5, 12, 6,  &
       6, 6, 6, 5, 6, 6, 6, 12, 8 ]

  real(rkx), dimension(22), parameter :: cdepuv = &
      [ 100.0_rkx, 100.0_rkx, 100.0_rkx, 100.0_rkx, 100.0_rkx, &
        100.0_rkx, 100.0_rkx, 100.0_rkx, 100.0_rkx, 100.0_rkx, &
        100.0_rkx, 100.0_rkx, 100.0_rkx, 100.0_rkx, 100.0_rkx, &
        100.0_rkx, 100.0_rkx, 100.0_rkx, 100.0_rkx, 100.0_rkx, &
        100.0_rkx, 100.0_rkx ]
  !
  ! rough is an aerodynamic roughness length (m) =approx 0.1*veg*height
  ! also used snow masking depth in subrout albedo

  ! Graziano - 2018-02-09 - switch to BATS defaults.
  real(rkx), dimension(22), parameter :: crough = &
    [ 0.1000_rkx, 0.0300_rkx, 1.0000_rkx, 1.0000_rkx, &
      1.0000_rkx, 1.0000_rkx, 0.3000_rkx, 0.0050_rkx, &
      0.0300_rkx, 0.1000_rkx, 0.0300_rkx, 0.0050_rkx, &
      0.1000_rkx, 0.0002_rkx, 0.0004_rkx, 0.2500_rkx, &
      0.1000_rkx, 1.0000_rkx, 0.5000_rkx, 0.3000_rkx, &
      2.0000_rkx, 1.0000_rkx ]

  ! real(rkx), dimension(22), parameter :: crough = &
  !  [ 0.0800_rkx, 0.0500_rkx, 1.0000_rkx, 1.0000_rkx, &
  !    0.8000_rkx, 2.0000_rkx, 0.1000_rkx, 0.0500_rkx, &
  !    0.0400_rkx, 0.0600_rkx, 0.1000_rkx, 0.0100_rkx, &
  !    0.0300_rkx, 0.0004_rkx, 0.0004_rkx, 0.1000_rkx, &
  !    0.0100_rkx, 0.8000_rkx, 0.3000_rkx, 0.3000_rkx, &
  !    1.5000_rkx, 0.4000_rkx ]

  ! Graziano - 2018-02-09 - Use values from
  !
  ! Leiming Zhang, Sunling Gong, Jacob Padro, Len Barrie
  ! A size-segregated particle dry deposition scheme for an atmospheric
  ! aerosol module, https://doi.org/10.1016/S1352-2310(00)00326-5
  !
  ! with the following lookup table from BATS classes into Table 3 from
  ! the paper:
  !
  ! BATS   :  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
  ! TABLE 3:  7  6  1  3  4  2  6  8  9  7  8 12 11 13 14 10 10  5  7 11 15 15
  !
   real(rkx), dimension(22) , parameter :: aest = &
    [  1.20_rkx,  1.20_rkx,  1.00_rkx,  1.10_rkx,  0.80_rkx, &
       0.60_rkx,  1.20_rkx, 50.00_rkx, 50.00_rkx,  1.20_rkx, &
      50.00_rkx, 50.00_rkx,  2.00_rkx ,100.00_rkx ,100.00_rkx, &
       1.30_rkx,  1.30_rkx,  0.80_rkx,  1.20_rkx,  2.00_rkx, &
       1.50_rkx,  1.50_rkx ]

  ! real(rkx), dimension(22) , parameter :: aest = &
  !  [  0.80_rkx,  0.80_rkx, 0.80_rkx, 0.80_rkx, 1.20_rkx, &
  !     1.20_rkx,  2.00_rkx, 1.50_rkx, 1.50_rkx, 2.00_rkx, &
  !    15.00_rkx, 15.00_rkx, 1.50_rkx, 1.50_rkx, 1.50_rkx, &
  !    15.00_rkx,  1.20_rkx, 1.20_rkx, 1.20_rkx, 1.20_rkx, &
  !     1.20_rkx,  1.20_rkx ]

   real(rkx), dimension(22) , parameter :: agam = &
     [ 0.5400_rkx,  0.5400_rkx,  0.5600_rkx,  0.5600_rkx, &
       0.5600_rkx,  0.5800_rkx,  0.5400_rkx,  0.5400_rkx, &
       0.5400_rkx,  0.5400_rkx,  0.5400_rkx,  0.5400_rkx, &
       0.5400_rkx,  0.5000_rkx,  0.5000_rkx,  0.5400_rkx, &
       0.5400_rkx,  0.5600_rkx,  0.5400_rkx,  0.5400_rkx, &
       0.5600_rkx,  0.5600_rkx ]
  !
  ! BATS   :  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
  ! TABLE 3:  7  6  1  3  4  2  6  8  9  7  8 12 11 13 14 10 10  5  7 11 15 15
  !
   real(rkx), dimension(22) , parameter :: arye = &
     [  2.0000_rkx,  2.0000_rkx,  2.0000_rkx,  2.0000_rkx, &
        5.0000_rkx,  5.0000_rkx,  5.0000_rkx,  0.0001_rkx, &
        0.0001_rkx,  2.0000_rkx,  0.0001_rkx,  0.0001_rkx, &
       10.0000_rkx,  0.0001_rkx,  0.0001_rkx, 10.0000_rkx, &
       10.0000_rkx,  5.0000_rkx,  2.0000_rkx, 10.0000_rkx, &
       10.0000_rkx, 10.000_rkx ]


  ! From the paper above, the A parameter used in Stokes number computation
  ! for vegetated classes.

   real(rkx), dimension(22) , parameter :: ast = &
      [  3.0_rkx,  3.0_rkx,  2.0_rkx,  3.0_rkx,  8.0_rkx,  &
         5.0_rkx,  2.0_rkx,  3.0_rkx,  0.0_rkx,  0.0_rkx,  &
        10.0_rkx,  0.0_rkx,  0.0_rkx, 10.0_rkx,  0.0_rkx,  &
         0.0_rkx, 10.0_rkx,  5.0_rkx,  3.0_rkx, 10.0_rkx,  &
        10.0_rkx, 10.0_rkx ]

end module mod_che_param
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
