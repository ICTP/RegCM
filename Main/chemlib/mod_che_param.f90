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

  public
!
  integer(ik4) :: ichaer

  real(rk8) , pointer , dimension(:) :: chlevs
!
  !FAB : temporaire car defini dans les modules de deposts
  real(rk8) , dimension(22) :: aest , arye
  ! Stokes parameters

  real(rk8) , dimension(12) :: cxmopor
  real(rk8) , dimension(22) :: cdepuv , crough
  integer(ik4) , dimension(22) :: ciexsol
  !
  !FAB : redefine soil prop. similar to bats for chemistry externalisation.
  ! think about an interface! 
  !
  !     ******      xmopor is fraction of soil that is voids
  data cxmopor/0.33D0 , 0.36D0 , 0.39D0 , 0.42D0 , 0.45D0 , 0.48D0 , &
               0.51D0 , 0.54D0 , 0.57D0 , 0.60D0 , 0.63D0 , 0.66D0/
  data ciexsol/6 , 6 , 6 , 6 , 7 , 8 , 6 , 1 , 6 , 6 , 5 , 12, 6 ,  &
               6 , 6 , 6 , 5 , 6 , 6 , 6 , 12, 8/
  data cdepuv /22*100.0D0/
  !
  ! rough is an aerodynamic roughness length (m) =approx 0.1*veg*height
  ! also used snow masking depth in subrout albedo
  data crough /0.08D0 , 0.05D0 , 2*1.0D0 , 0.8D0 , 2.0D0  , 0.1D0  , &
               0.05D0 , 0.04D0 , 0.06D0 ,  0.1D0 , 0.01D0 , 0.03D0 , &
               2*0.0004D0 , 2*0.1D0 , 0.8D0 , 2*0.3D0, 1.5D0, 0.40D0 /

  data aest /0.80D0 , 0.80D0 , 0.8D0 , 0.8D0 , 1.2D0 , 1.20D0 , &
       2.0D0 , 1.5D0 ,  1.5D0 , 2.0D0 , 15.0D0 , 15.0D0 , 1.5D0 ,   &
       1.5D0 , 1.5D0 , 15.0D0 , 1.2D0 , 1.2D0 , 1.2D0 , 1.2D0 ,     &
       1.2D0 , 1.2D0 /
!
  data arye /0.5D0 , 5.0D0 , 0.5D0 , 5.0D0 , 1.0D0 , 1.0D0 ,    &
     0.0001D0 , 5.0D0 , 10.0D0 , 10.0D0 , 0.0001D0 , 0.0001D0 ,     &
     0.56D0 , 0.56D0 , 0.56D0 , 0.56D0 ,  0.56D0 , 0.56D0 , 0.56D0 ,&
     0.56D0 , 1.0D0 , 1.0D0 /
!
end module mod_che_param
