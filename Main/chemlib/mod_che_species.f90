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

    integer , parameter :: ninchcodes = 46
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

    real , parameter :: w_no2 = 46.0D0
    real , parameter :: w_no  = 30.0D0

    real , parameter :: w_hono = 47.0D0
    real , parameter :: w_no3  = 62.0D0
    real , parameter :: w_n2o5 = 108.0D0
    real , parameter :: w_hno4 = 79.0D0
    real , parameter :: w_hno3 = 63.0D0
    real , parameter :: w_o3   = 48.0D0
    real , parameter :: w_h2o2 = 34.0D0
    real , parameter :: w_h2o  = 18.0D0

    real , parameter :: w_so2  = 64.0D0
    real , parameter :: w_sulf = 98.0D0
    real , parameter :: w_co   = 28.0D0
    real , parameter :: w_co2  = 44.0D0
    real , parameter :: w_h2   = 2.0D0

    real , parameter :: w_oh   = 17.0D0
    real , parameter :: w_ho2  = 33.0D0

    real , parameter :: w_ch4    = 16.0D0
    real , parameter :: w_ethan  = 30.0D0
    real , parameter :: w_hc3    = 44.0D0
    real , parameter :: w_hc5    = 72.0D0
    real , parameter :: w_hc8    = 114.0D0

    real , parameter :: w_ethene = 28.0D0
    real , parameter :: w_ol2    = 28.0D0
    real , parameter :: w_olt    = 42.0D0
    real , parameter :: w_prpe   = 42.0D0
    real , parameter :: w_bute   = 56.0D0
    real , parameter :: w_oli    = 56.0D0
    real , parameter :: w_isop   = 68.0D0

    real , parameter :: w_tolu   = 92.0D0
    real , parameter :: w_csl    = 108.0D0
    real , parameter :: w_xyle   = 106.0D0

    real , parameter :: w_hcho    = 30.0D0
    real , parameter :: w_ald2    = 44.0D0
    real , parameter :: w_ket     = 72.0D0
    real , parameter :: w_onit    = 119.0D0
    real , parameter :: w_gly     = 58.0D0
    real , parameter :: w_mgly    = 72.0D0

    real , parameter :: w_pan     = 121.0D0
    real , parameter :: w_tpan    = 147.0D0
    real , parameter :: w_alco    = 32.0D0

    real , parameter :: w_dms     = 62.0D0
    real , parameter :: w_rooh    = 48.0D0
    real , parameter :: w_nh3     = 17.0D0
    real , parameter :: w_c2h6    = 30.07D0
    real , parameter :: w_c3h8    = 44.1D0
    real , parameter :: w_alk4    = 58.12D0
    real , parameter :: w_alk7    = 100.20D0
    real , parameter :: w_mo2     = 47.0D0
    real , parameter :: w_acet    = 58.08D0
    real , parameter :: w_moh     = 32.04D0
    real , parameter :: w_eoh     = 46.07D0
    real , parameter :: w_benz    = 78.11D0

  contains

    subroutine allocate_che_species_bc
      call getmem4d(chbc0%bc,1,iy,1,kz,1,jxp,1,nchsp,'mod_che_species:chbc0')
      call getmem4d(chbc1%bc,1,iy,1,kz,1,jxp,1,nchsp,'mod_che_species:chbc1')
      call getmem4d(oxbc0%bc,1,iy,1,kz,1,jxp,1,noxsp,'mod_che_species:oxbc0')
    end subroutine allocate_che_species_bc

  end module mod_che_species
