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

module mod_isaatm
  use mod_realkinds
  use mod_constants

  implicit none (type, external)
  private

  public :: isaatm

  contains

  subroutine isaatm(height,temperature,pressure,density)
    implicit none (type, external)
    real(rkx), intent(in) :: height       ! Height in meter
    real(rkx), intent(out) :: temperature ! Temperature in K
    real(rkx), intent(out) :: pressure    ! Pressure in Pa
    real(rkx), intent(out) :: density     ! Density in kg/m^3

    real(rkx) :: h, i, m, n, o

    real(rkx), parameter :: height0 = 0.0_rkx
    real(rkx), parameter :: height1 = 11000.0_rkx
    real(rkx), parameter :: height2 = 20000.0_rkx
    real(rkx), parameter :: height3 = 32000.0_rkx
    real(rkx), parameter :: height4 = 47000.0_rkx

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
