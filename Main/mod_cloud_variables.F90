!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by
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

module mod_cloud_variables
   
    use mod_realkinds
    use mod_constants



    implicit none
     
    !numbers
    real(rk8) , parameter :: rcldiff = 3.D-6
    logical , parameter :: lmicro = .true.
    logical , parameter :: fscheme = .true.
    real(rk8) , parameter :: ksemi =  d_one
    real(rk8) , parameter :: rlcritsnow = 3.D-5   !critical autoconversion
    !real(rk8) , parameter :: zauto_rate_khair = 0.355D0 ! microphysical terms
    real(rk8) , parameter :: zauto_expon_khair = 1.47D0
    real(rk8) , parameter :: zrldcp = d_one/(wlhsocp-wlhvocp)  ! Cp/Lf
   ! 1/autoconversion time scale (s)
    !real(rk8) , parameter :: rkconv = d_one/6000.0D0
    ! real(rk8) , parameter :: zauto_rate_sundq = 0.5D-3
    !real(rk8) , parameter :: zauto_rate_kessl = 1.D-3    !giusto!
    real(rk8) , parameter :: zautocrit_kessl = 5.D-4          !giusto!
    !real(rk8) , parameter :: zauto_rate_klepi = 0.5D-3
    real(rk8) , parameter :: zepsec = 1.D-14
    ! real(rk8) , parameter :: qs0 = 1.0D-3             !g g^(-1)
    !real(rk8) , parameter :: rcovpmin = 0.1
    real(rk8) , parameter :: rclcrit_land = 5.D-4
    real(rk8) , parameter :: rclcrit_sea = 3.D-4
    real(rk8) , parameter :: rprc1 = 3.D2           ! in Sundqvist = 300
    real(rk8) , parameter :: ramid = 0.8D0
    !real(rk8) , parameter :: rlmin = 1.D-8
    ! max threshold rh for evaporation
    ! for a precip coverage of zero rmfsens20081111
    real(rk8) , parameter :: rprecrhmax = 0.7D0
    !evaporation rate coefficient
    ! Numerical fit to wet bulb temperature
    real(rk8) , parameter :: ztw1 = 1329.31D0
    real(rk8) , parameter :: ztw2 = 0.0074615D0
    real(rk8) , parameter :: ztw3 = 0.85D5
    real(rk8) , parameter :: ztw4 = 40.637D0
    real(rk8) , parameter :: ztw5 = 275.0D0
    real(rk8) , parameter :: rtaumel = 1.1880D4
    ! variables/constants for the supersaturation
    real(rk8) , parameter :: r5les =  4097.9337D0  !r3les*(rtt-r4les)
    real(rk8) , parameter :: r5ies =  5807.547D0 !r3ies*(rtt-r4ies)
   real(rk8) , parameter :: r4les =  35.86D0
    real(rk8) , parameter :: r4ies = 7.66D0
    ! temperature homogeneous freezing
    real(rk8) , parameter :: thomo = 235.16D0  ! -38.00 Celsius
    ! Cloud fraction threshold that defines cloud top
    real(rk8), parameter :: rcldtopcf = d_r100
    ! Fraction of deposition rate in cloud top layer
    real(rk8), parameter :: rdepliqrefrate = d_r10
    ! Depth of supercooled liquid water layer (m)
    real(rk8), parameter :: rdepliqrefdepth = 500.0D0
    ! initial mass of ice particle
    real(rk8), parameter :: riceinit = 1.D-12

end module mod_cloud_variables
