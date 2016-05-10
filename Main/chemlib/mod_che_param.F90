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
  real(rkx) , dimension(22) :: aest , arye
  ! Stokes parameters

  real(rkx) , dimension(12) :: cxmopor
  real(rkx) , dimension(22) :: cdepuv , crough
  integer(ik4) , dimension(22) :: ciexsol
  !
  !FAB : redefine soil prop. similar to bats for chemistry externalisation.
  ! think about an interface!
  !
  !     ******      xmopor is fraction of soil that is voids
  data cxmopor/0.33_rkx , 0.36_rkx , 0.39_rkx , 0.42_rkx , 0.45_rkx , 0.48_rkx , &
               0.51_rkx , 0.54_rkx , 0.57_rkx , 0.60_rkx , 0.63_rkx , 0.66_rkx/
  data ciexsol/6 , 6 , 6 , 6 , 7 , 8 , 6 , 1 , 6 , 6 , 5 , 12, 6 ,  &
               6 , 6 , 6 , 5 , 6 , 6 , 6 , 12, 8/
  data cdepuv /22*100.0_rkx/
  !
  ! rough is an aerodynamic roughness length (m) =approx 0.1*veg*height
  ! also used snow masking depth in subrout albedo
  data crough /0.08_rkx , 0.05_rkx , 2*1.0_rkx , 0.8_rkx , 2.0_rkx  , 0.1_rkx  , &
               0.05_rkx , 0.04_rkx , 0.06_rkx ,  0.1_rkx , 0.01_rkx , 0.03_rkx , &
               2*0.0004_rkx , 2*0.1_rkx , 0.8_rkx , 2*0.3_rkx, 1.5_rkx, 0.40_rkx /

  data aest /0.80_rkx , 0.80_rkx , 0.8_rkx , 0.8_rkx , 1.2_rkx , 1.20_rkx , &
       2.0_rkx , 1.5_rkx ,  1.5_rkx , 2.0_rkx , 15.0_rkx , 15.0_rkx , 1.5_rkx ,   &
       1.5_rkx , 1.5_rkx , 15.0_rkx , 1.2_rkx , 1.2_rkx , 1.2_rkx , 1.2_rkx ,     &
       1.2_rkx , 1.2_rkx /
!
  data arye /0.5_rkx , 5.0_rkx , 0.5_rkx , 5.0_rkx , 1.0_rkx , 1.0_rkx ,    &
     0.0001_rkx , 5.0_rkx , 10.0_rkx , 10.0_rkx , 0.0001_rkx , 0.0001_rkx ,     &
     0.56_rkx , 0.56_rkx , 0.56_rkx , 0.56_rkx ,  0.56_rkx , 0.56_rkx , 0.56_rkx ,&
     0.56_rkx , 1.0_rkx , 1.0_rkx /
!
end module mod_che_param
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
