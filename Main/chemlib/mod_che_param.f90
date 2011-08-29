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

  use m_realkinds

  public
!
  logical :: lch
  integer :: ichdir

  real(dp) :: chfrq , rafrq , mdfrq

  real(dp) , pointer , dimension(:) :: chlevs
!
  real(8) , dimension(22) :: aest , arye
!
! Stokes parameters
!
  data aest     /0.80D0 , 0.80D0 , 0.8D0 , 0.8D0 , 1.2D0 , 1.20D0 , &
       2.0D0 , 1.5D0 ,  1.5D0 , 2.0D0 , 15.0D0 , 15.0D0 , 1.5D0 ,   &
       1.5D0 , 1.5D0 , 15.0D0 , 1.2D0 , 1.2D0 , 1.2D0 , 1.2D0 ,     &
       1.2D0 , 1.2D0 /
!
  data arye     /0.5D0 , 5.0D0 , 0.5D0 , 5.0D0 , 1.0D0 , 1.0D0 ,    &
     0.0001D0 , 5.0D0 , 10.0D0 , 10.0D0 , 0.0001D0 , 0.0001D0 ,     &
     0.56D0 , 0.56D0 , 0.56D0 , 0.56D0 ,  0.56D0 , 0.56D0 , 0.56D0 ,&
     0.56D0 , 1.0D0 , 1.0D0 /
!
end module mod_che_param
