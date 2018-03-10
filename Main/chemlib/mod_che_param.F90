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

module mod_che_param

  use mod_intkinds
  use mod_realkinds

  implicit none

  public

  integer(ik4) :: ichaer

  !FAB : temporaire car defini dans les modules de deposts
  real(rkx) , dimension(22) :: aest , arye , ast , agam
  ! Stokes parameters

  real(rkx) , dimension(12) :: cxmopor
  real(rkx) , dimension(22) :: cdepuv , crough
  integer(ik4) , dimension(22) :: ciexsol
  !
  !FAB : redefine soil prop. similar to bats for chemistry externalisation.
  ! think about an interface!
  !
  !     ******      xmopor is fraction of soil that is voids
  data cxmopor / 0.33_rkx , 0.36_rkx , 0.39_rkx , 0.42_rkx , 0.45_rkx , &
                 0.48_rkx , 0.51_rkx , 0.54_rkx , 0.57_rkx , 0.60_rkx , &
                 0.63_rkx , 0.66_rkx/

  data ciexsol/6 , 6 , 6 , 6 , 7 , 8 , 6 , 1 , 6 , 6 , 5 , 12, 6 ,  &
               6 , 6 , 6 , 5 , 6 , 6 , 6 , 12, 8/

  data cdepuv /22*100.0_rkx/
  !
  ! rough is an aerodynamic roughness length (m) =approx 0.1*veg*height
  ! also used snow masking depth in subrout albedo

  ! Graziano - 2018-02-09 - switch to BATS defaults.
  data crough / 0.1000_rkx , 0.0300_rkx , 1.0000_rkx , 1.0000_rkx , &
                1.0000_rkx , 1.0000_rkx , 0.3000_rkx , 0.0050_rkx , &
                0.0300_rkx , 0.1000_rkx , 0.0300_rkx , 0.0050_rkx , &
                0.1000_rkx , 0.0002_rkx , 0.0004_rkx , 0.2500_rkx , &
                0.1000_rkx , 1.0000_rkx , 0.5000_rkx , 0.3000_rkx , &
                2.0000_rkx , 1.0000_rkx /

  !data crough /0.0800_rkx , 0.0500_rkx , 1.0000_rkx , 1.0000_rkx , &
  !             0.8000_rkx , 2.0000_rkx , 0.1000_rkx , 0.0500_rkx , &
  !             0.0400_rkx , 0.0600_rkx , 0.1000_rkx , 0.0100_rkx , &
  !             0.0300_rkx , 0.0004_rkx , 0.0004_rkx , 0.1000_rkx , &
  !             0.0100_rkx , 0.8000_rkx , 0.3000_rkx , 0.3000_rkx , &
  !             1.5000_rkx , 0.4000_rkx /

  ! Graziano - 2018-02-09 - Use values from
  !
  ! Leiming Zhang, Sunling Gong, Jacob Padro, Len Barrie
  ! A size-segregated particle dry deposition scheme for an atmospheric
  ! aerosol module , https://doi.org/10.1016/S1352-2310(00)00326-5
  !
  ! with the following lookup table from BATS classes into Table 3 from
  ! the paper:
  !
  ! BATS   :  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
  ! TABLE 3:  7  6  1  3  4  2  6  8  9  7  8 12 11 13 14 10 10  5  7 11 15 15
  !
  data aest  /  1.20_rkx ,  1.20_rkx ,  1.00_rkx ,  1.10_rkx ,  0.80_rkx , &
                0.60_rkx ,  1.20_rkx , 50.00_rkx , 50.00_rkx ,  1.20_rkx , &
               50.00_rkx , 50.00_rkx ,  2.00_rkx ,100.00_rkx ,100.00_rkx , &
                1.30_rkx ,  1.30_rkx ,  0.80_rkx ,  1.20_rkx ,  2.00_rkx , &
                1.50_rkx ,  1.50_rkx /

  !data aest /  0.80_rkx ,  0.80_rkx , 0.80_rkx , 0.80_rkx , 1.20_rkx , &
  !             1.20_rkx ,  2.00_rkx , 1.50_rkx , 1.50_rkx , 2.00_rkx , &
  !            15.00_rkx , 15.00_rkx , 1.50_rkx , 1.50_rkx , 1.50_rkx , &
  !            15.00_rkx ,  1.20_rkx , 1.20_rkx , 1.20_rkx , 1.20_rkx , &
  !             1.20_rkx ,  1.20_rkx /

  data agam /  0.5400_rkx ,  0.5400_rkx ,  0.5600_rkx ,  0.5600_rkx , &
               0.5600_rkx ,  0.5800_rkx ,  0.5400_rkx ,  0.5400_rkx , &
               0.5400_rkx ,  0.5400_rkx ,  0.5400_rkx ,  0.5400_rkx , &
               0.5400_rkx ,  0.5000_rkx ,  0.5000_rkx ,  0.5400_rkx , &
               0.5400_rkx ,  0.5600_rkx ,  0.5400_rkx ,  0.5400_rkx , &
               0.5600_rkx ,  0.5600_rkx /

  !
  ! BATS   :  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
  ! TABLE 3:  7  6  1  3  4  2  6  8  9  7  8 12 11 13 14 10 10  5  7 11 15 15
  !
   data arye / 2.0000_rkx ,  2.0000_rkx ,  2.0000_rkx ,  2.0000_rkx , &
               5.0000_rkx ,  5.0000_rkx ,  5.0000_rkx ,  0.0001_rkx , &
               0.0001_rkx ,  2.0000_rkx ,  0.0001_rkx ,  0.0001_rkx , &
              10.0000_rkx ,  0.0001_rkx ,  0.0001_rkx , 10.0000_rkx , &
              10.0000_rkx ,  5.0000_rkx ,  2.0000_rkx , 10.0000_rkx , &
              10.0000_rkx , 10.000_rkx /


  ! From the paper above, the A parameter used in Stokes number computation
  ! for vegetated classes.

  data ast /  3.0 ,  3.0 ,  2.0 ,  3.0 ,  8.0 ,  5.0 ,  2.0 ,  3.0 , &
              0.0 ,  0.0 , 10.0 ,  0.0 ,  0.0 , 10.0 ,  0.0 ,  0.0 , &
             10.0 ,  5.0 ,  3.0 , 10.0 , 10.0 , 10.0 /

end module mod_che_param
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
