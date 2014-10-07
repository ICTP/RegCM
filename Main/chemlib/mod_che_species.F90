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

module mod_che_species

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil

  implicit none

  public

  integer(ik4) , parameter :: nchsp = 25
  integer(ik4) , parameter :: noxsp = 5

  integer(ik4) , parameter :: ninchcodes = 46
  character(len=8) , dimension(ninchcodes) :: inchtrnames =  &
   (/'O3      ','NO      ','NO2     ','NO3     ','OH      ', &
     'HO2     ','H2O2    ','HNO2    ','HNO3    ','HNO4    ', &
     'SULF    ','SO4     ','H2SO4   ','HONO    ','N2O5    ', &
     'HC      ','HCR     ','C2H4    ','OLT     ','OLI     ', &
     'ALK4    ','ALK7    ','CO      ','HCHO    ','ALD2    ', &
     'ETHENE  ','C2H6    ','C3H8    ','ISOP    ','TOLUENE ', &
     'XYL     ','NH3     ','PAN     ','ROOH    ','ACET    ', &
     'BENZ    ','CH4     ','MOH     ','EOH     ','ACO2    ', &
     'CO2     ','DMS     ','NOX     ','HOX     ','SOX     ', &
     'PAR     '/)

  integer(ik4) , parameter :: ich_o3       = 1
  integer(ik4) , parameter :: ich_no       = 2
  integer(ik4) , parameter :: ich_no2      = 3
  integer(ik4) , parameter :: ich_hno3     = 4
  integer(ik4) , parameter :: ich_n2o5     = 5
  integer(ik4) , parameter :: ich_h2o2     = 6
  integer(ik4) , parameter :: ich_ch4      = 7
  integer(ik4) , parameter :: ich_co       = 8
  integer(ik4) , parameter :: ich_hcho     = 9
  integer(ik4) , parameter :: ich_ch3oh    = 10
  integer(ik4) , parameter :: ich_c2h5oh   = 11
  integer(ik4) , parameter :: ich_c2h4     = 12
  integer(ik4) , parameter :: ich_c2h6     = 13
  integer(ik4) , parameter :: ich_ch3cho   = 14
  integer(ik4) , parameter :: ich_c3h6     = 15
  integer(ik4) , parameter :: ich_c3h8     = 16
  integer(ik4) , parameter :: ich_ch3coch3 = 17
  integer(ik4) , parameter :: ich_bigene   = 18
  integer(ik4) , parameter :: ich_bigalk   = 19
  integer(ik4) , parameter :: ich_isop     = 20
  integer(ik4) , parameter :: ich_tolue    = 21
  integer(ik4) , parameter :: ich_pan      = 22
  integer(ik4) , parameter :: ich_so2      = 23
  integer(ik4) , parameter :: ich_dms      = 24
  integer(ik4) , parameter :: ich_so4      = 25

  integer(ik4) , parameter :: iox_oh       = 1
  integer(ik4) , parameter :: iox_ho2      = 2
  integer(ik4) , parameter :: iox_o3       = 3
  integer(ik4) , parameter :: iox_no3      = 4
  integer(ik4) , parameter :: iox_h2o2     = 5

end module mod_che_species
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
