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

    use m_realkinds
    use mod_dynparam
    use mod_memutil

    public

    integer , parameter :: nchsp = 25
    integer , parameter :: noxsp = 5

    integer , parameter :: ich_o3       = 1
    integer , parameter :: ich_no       = 2
    integer , parameter :: ich_no2      = 3
    integer , parameter :: ich_hno3     = 4
    integer , parameter :: ich_n2o5     = 5
    integer , parameter :: ich_h2o2     = 6
    integer , parameter :: ich_ch4      = 7
    integer , parameter :: ich_co       = 8
    integer , parameter :: ich_hcho     = 9
    integer , parameter :: ich_ch3oh    = 10
    integer , parameter :: ich_c2h5oh   = 11
    integer , parameter :: ich_c2h4     = 12
    integer , parameter :: ich_c2h6     = 13
    integer , parameter :: ich_ch3cho   = 14
    integer , parameter :: ich_c3h6     = 15
    integer , parameter :: ich_c3h8     = 16
    integer , parameter :: ich_ch3coch3 = 17
    integer , parameter :: ich_bigene   = 18
    integer , parameter :: ich_bigalk   = 19
    integer , parameter :: ich_isop     = 20
    integer , parameter :: ich_tolue    = 21
    integer , parameter :: ich_pan      = 22
    integer , parameter :: ich_so2      = 23
    integer , parameter :: ich_dms      = 24
    integer , parameter :: ich_so4      = 25

    integer , parameter :: iox_oh       = 1
    integer , parameter :: iox_ho2      = 2
    integer , parameter :: iox_o3       = 3
    integer , parameter :: iox_no3      = 4
    integer , parameter :: iox_h2o2     = 5

    type che_sp_bc
      real(dp) , pointer , dimension(:,:,:,:) :: bc
    end type

    type che_ox_bc
      real(dp) , pointer , dimension(:,:,:,:) :: bc
    end type

    type(che_sp_bc) :: chbc0 , chbc1
    type(che_ox_bc) :: oxbc0

  contains

    subroutine allocate_che_species_bc
      call getmem4d(chbc0%bc,1,iy,1,kz,1,jxp,1,nchsp,'mod_che_species:chbc0')
      call getmem4d(chbc1%bc,1,iy,1,kz,1,jxp,1,nchsp,'mod_che_species:chbc1')
      call getmem4d(oxbc0%bc,1,iy,1,kz,1,jxp,1,noxsp,'mod_che_species:oxbc0')
    end subroutine allocate_che_species_bc

  end module mod_che_species
