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

module mod_cu_tables

  use m_realkinds
  use mod_constants

  implicit none

  private

  public :: jptlucu1        ! lookup table lower bound
  public :: jptlucu2        ! lookup table upper bound
  public :: tlucua          ! table -- e_s*rgas/rwat
  public :: tlucub          ! table -- for derivative calculation: d es/ d t
  public :: tlucuc          ! table -- l/cp
  public :: tlucuaw         ! table
  public :: lookupoverflow  ! lookup table overflow flag

  ! subroutines public

  public :: init_convect_tables ! initialize luts 

  integer, parameter :: jptlucu1 =  50000  ! lookup table lower bound
  integer, parameter :: jptlucu2 = 400000  ! lookup table upper bound

  logical :: lookupoverflow = .false.      ! preset with false
  
  real(dp) :: tlucua(jptlucu1:jptlucu2)    ! table - e_s*rgas/rwat
  real(dp) :: tlucub(jptlucu1:jptlucu2)    ! table - for derivative calculation
  real(dp) :: tlucuc(jptlucu1:jptlucu2)    ! table - l/cp
  real(dp) :: tlucuaw(jptlucu1:jptlucu2)   ! table

!------------------------------------------------------------------------------

contains

  subroutine init_convect_tables

    real(dp), parameter :: zavl1 = -6096.9385_dp
    real(dp), parameter :: zavl2 =    21.2409642_dp
    real(dp), parameter :: zavl3 =    -2.711193_dp
    real(dp), parameter :: zavl4 =     1.673952_dp
    real(dp), parameter :: zavl5 =     2.433502_dp 

    real(dp), parameter :: zavi1 = -6024.5282_dp
    real(dp), parameter :: zavi2 =    29.32707_dp
    real(dp), parameter :: zavi3 =     1.0613868_dp
    real(dp), parameter :: zavi4 =    -1.3198825_dp
    real(dp), parameter :: zavi5 =    -0.49382577_dp        

    real(dp) :: z5alvcp, z5alscp, zalvdcp, zalsdcp
    real(dp) :: ztt, zldcp
    real(dp) :: zcvm3, zcvm4, zcvm5
    real(dp) :: zavm1, zavm2, zavm3, zavm4, zavm5

    integer :: it

    z5alvcp = c5les*wlhv/cpd
    z5alscp = c5ies*wlhs/cpd

    zalvdcp = wlhv/cpd
    zalsdcp = wlhs/cpd

    do it = jptlucu1, jptlucu2
      ztt = 0.001_dp*it
      if ((ztt-tzero) > 0.0_dp) then
        zcvm3 = c3les
        zcvm4 = c4les
        zcvm5 = z5alvcp
        zldcp = zalvdcp
        zavm1 = zavl1
        zavm2 = zavl2
        zavm3 = zavl3
        zavm4 = zavl4
        zavm5 = zavl5
      else
        zcvm3 = c3ies
        zcvm4 = c4ies
        zcvm5 = z5alscp
        zldcp = zalsdcp
        zavm1 = zavi1
        zavm2 = zavi2
        zavm3 = zavi3
        zavm4 = zavi4
        zavm5 = zavi5
      end if
      tlucuc(it)  = zldcp
      tlucua(it)  = exp((zavm1/ztt+zavm2+zavm3*0.01_dp* &
                    ztt+zavm4*ztt*ztt*1.e-5_dp+zavm5*log(ztt)))*rgas/rwat
      tlucub(it)  = zcvm5*(1.0_dp/(ztt-zcvm4))**2
      tlucuaw(it) = exp((zavl1/ztt+zavl2+zavl3*0.01_dp* &
                    ztt+zavl4*ztt*ztt*1.e-5_dp+zavl5*log(ztt)))*rgas/rwat
    end do
    
  end subroutine init_convect_tables

end module mod_cu_tables
