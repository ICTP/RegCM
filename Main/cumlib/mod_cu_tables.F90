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

  use mod_intkinds
  use mod_realkinds
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

  real(rkx) :: tlucua(jptlucu1:jptlucu2)    ! table - e_s*rgas/rwat
  real(rkx) :: tlucub(jptlucu1:jptlucu2)    ! table - for derivative calculation
  real(rkx) :: tlucuc(jptlucu1:jptlucu2)    ! table - l/cp
  real(rkx) :: tlucuaw(jptlucu1:jptlucu2)   ! table

  contains

  subroutine init_convect_tables
    implicit none
    real(rkx), parameter :: zavl1 = -6096.9385_rkx
    real(rkx), parameter :: zavl2 =    21.2409642_rkx
    real(rkx), parameter :: zavl3 =    -2.711193_rkx
    real(rkx), parameter :: zavl4 =     1.673952_rkx
    real(rkx), parameter :: zavl5 =     2.433502_rkx

    real(rkx), parameter :: zavi1 = -6024.5282_rkx
    real(rkx), parameter :: zavi2 =    29.32707_rkx
    real(rkx), parameter :: zavi3 =     1.0613868_rkx
    real(rkx), parameter :: zavi4 =    -1.3198825_rkx
    real(rkx), parameter :: zavi5 =    -0.49382577_rkx

    real(rkx) :: z5alvcp, z5alscp, zalvdcp, zalsdcp
    real(rkx) :: ztt, zldcp
    real(rkx) :: zcvm3, zcvm4, zcvm5
    real(rkx) :: zavm1, zavm2, zavm3, zavm4, zavm5

    integer(ik4) :: it

    z5alvcp = c5les*wlhv*rcpd
    z5alscp = c5ies*wlhs*rcpd

    zalvdcp = wlhv*rcpd
    zalsdcp = wlhs*rcpd

    do it = jptlucu1, jptlucu2
      ztt = 0.001_rkx*it
      if ((ztt-tzero) > 0.0_rkx) then
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
      tlucua(it)  = exp((zavm1/ztt+zavm2+zavm3*0.01_rkx* &
                    ztt+zavm4*ztt*ztt*1.e-5_rkx+zavm5*log(ztt)))*rgas/rwat
      tlucub(it)  = zcvm5*(1.0_rkx/(ztt-zcvm4))**2.0_rkx
      tlucuaw(it) = exp((zavl1/ztt+zavl2+zavl3*0.01_rkx* &
                    ztt+zavl4*ztt*ztt*1.e-5_rkx+zavl5*log(ztt)))*rgas/rwat
    end do

  end subroutine init_convect_tables

end module mod_cu_tables
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
