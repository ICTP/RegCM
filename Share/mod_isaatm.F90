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

module mod_isaatm
  use mod_realkinds
  use mod_constants

  implicit none
  private

  public :: isaatm

  contains

  subroutine isaatm(height,temperature,pressure,density)
    implicit none
    real(rkx) , intent(in) :: height       ! Height in meter
    real(rkx) , intent(out) :: temperature ! Temperature in K
    real(rkx) , intent(out) :: pressure    ! Pressure in Pa
    real(rkx) , intent(out) :: density     ! Density in kg/m^3

    real(rkx) :: h , i , m , n , o

    real(rkx) , parameter :: height0 = 0.0_rkx
    real(rkx) , parameter :: height1 = 11000.0_rkx
    real(rkx) , parameter :: height2 = 20000.0_rkx
    real(rkx) , parameter :: height3 = 32000.0_rkx
    real(rkx) , parameter :: height4 = 47000.0_rkx

    if ( height >= height0 .and. height <= height1 ) then
      temperature = stdt - lrate*height
      pressure = stdp * (temperature/stdt)**(egrav/(rgas*lrate))
      density = pressure/(rgas*temperature)
    else if ( height > height1 .and. height <= height2 ) then
      temperature = stdt - lrate*height1
      i = stdp * (temperature/stdt)**(egrav/(rgas*lrate))
      pressure = i * exp(-egrav*(height-height1)/(rgas*temperature))
      density = pressure/(rgas*temperature)
    else if ( height > height2 .and. height <= height3 ) then
      h = stdt - lrate*height1
      i = stdp * (h/stdt)**(egrav/(rgas*lrate))
      m = i * exp((-egrav*(height2-height1)/(rgas*h)))
      temperature = h + 0.001_rkx*(height -height2)
      pressure = m * (temperature/h)**(-(egrav/(rgas*0.001)))
      density = pressure/(rgas*temperature)
    else if ( height > height3 .and. height <= height4 ) then
      h = stdt - lrate*height1
      i = stdp * (h/stdt)**((egrav/(rgas*lrate)))
      m = i * exp(-egrav*(height2-height1)/(rgas*h))
      n = h + 0.001_rkx*(height3-height2)
      o = m * (n/h)**(-(egrav/(rgas*0.001_rkx)))
      n = h + 0.001_rkx*(height3-height2)
      temperature = n + 0.0028_rkx*(height-height3)
      pressure = o * (height/n)**(-(egrav/(rgas*0.0028_rkx)))
      density = pressure/(rgas/temperature)
    end if
  end subroutine isaatm

end module mod_isaatm
