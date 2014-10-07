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

module mod_ch_param

  use mod_intkinds
  use mod_realkinds

  public




!Mozart indcies
 INTEGER(ik4), PARAMETER :: mz_NO      = 1
 INTEGER(ik4), PARAMETER :: mz_NO2     = 2
 INTEGER(ik4), PARAMETER :: mz_N2O5    = 3
 INTEGER(ik4), PARAMETER :: mz_HNO3    = 4
 INTEGER(ik4), PARAMETER :: mz_HO2NO2  = 5
 INTEGER(ik4), PARAMETER :: mz_O3      = 6
 INTEGER(ik4), PARAMETER :: mz_H2O2    = 7
 INTEGER(ik4), PARAMETER :: mz_SO2     = 8
 INTEGER(ik4), PARAMETER :: mz_SO4     = 9
 INTEGER(ik4), PARAMETER :: mz_CH4     = 10
 INTEGER(ik4), PARAMETER :: mz_CH2O    = 11
 INTEGER(ik4), PARAMETER :: mz_CH3OH   = 12
 INTEGER(ik4), PARAMETER :: mz_PAN     = 13
 INTEGER(ik4), PARAMETER :: mz_C2H6    = 14
 INTEGER(ik4), PARAMETER :: mz_C3H8    = 15
 INTEGER(ik4), PARAMETER :: mz_BIGALK  = 16
 INTEGER(ik4), PARAMETER :: mz_C2H4    = 17
 INTEGER(ik4), PARAMETER :: mz_C3H6    = 18
 INTEGER(ik4), PARAMETER :: mz_BIGENE  = 19
 INTEGER(ik4), PARAMETER :: mz_TOLUENE = 20
 INTEGER(ik4), PARAMETER :: mz_ISOP    = 21
 INTEGER(ik4), PARAMETER :: mz_CH3CHO  = 22
 INTEGER(ik4), PARAMETER :: mz_CH3COOH = 23
 INTEGER(ik4), PARAMETER :: mz_GLYALD  = 24
 INTEGER(ik4), PARAMETER :: mz_CH3OOH  = 25
 INTEGER(ik4), PARAMETER :: mz_C2H5OOH = 26
 INTEGER(ik4), PARAMETER :: mz_CH3COCH3= 27
 INTEGER(ik4), PARAMETER :: mz_HYAC    = 28
 INTEGER(ik4), PARAMETER :: mz_CH3COCHO =29
 INTEGER(ik4), PARAMETER :: mz_ONIT    = 30
 INTEGER(ik4), PARAMETER :: mz_MEK     = 31
 INTEGER(ik4), PARAMETER :: mz_MVK     = 32
 INTEGER(ik4), PARAMETER :: mz_MACR    = 33
 INTEGER(ik4), PARAMETER :: mz_HYDRALD = 34
 INTEGER(ik4), PARAMETER :: mz_BIGALD  = 35
 INTEGER(ik4), PARAMETER :: mz_ISOPNO3 = 36
 INTEGER(ik4), PARAMETER :: mz_ONITR   = 37
 INTEGER(ik4), PARAMETER :: mz_CRESOL  = 38
 INTEGER(ik4), PARAMETER :: mz_CO      = 39
 INTEGER(ik4), PARAMETER :: mz_DMS     = 40


!CBMZ indcies
 INTEGER(ik4), PARAMETER :: cb_O3       = 1
 INTEGER(ik4), PARAMETER :: cb_NO       = 2
 INTEGER(ik4), PARAMETER :: cb_NO2      = 3
 INTEGER(ik4), PARAMETER :: cb_HNO3     = 4
 INTEGER(ik4), PARAMETER :: cb_HNO4     = 5
 INTEGER(ik4), PARAMETER :: cb_N2O5     = 6
 INTEGER(ik4), PARAMETER :: cb_H2O2     = 7
 INTEGER(ik4), PARAMETER :: cb_CH4      = 8
 INTEGER(ik4), PARAMETER :: cb_CO       = 9
 INTEGER(ik4), PARAMETER :: cb_SO2      = 10
 INTEGER(ik4), PARAMETER :: cb_H2SO4    = 11
 INTEGER(ik4), PARAMETER :: cb_DMS      = 12
 INTEGER(ik4), PARAMETER :: cb_PAR      = 13
 INTEGER(ik4), PARAMETER :: cb_C2H6     = 14
 INTEGER(ik4), PARAMETER :: cb_ETH      = 15
 INTEGER(ik4), PARAMETER :: cb_OLET     = 16
 INTEGER(ik4), PARAMETER :: cb_OLEI     = 17
 INTEGER(ik4), PARAMETER :: cb_TOL      = 18
 INTEGER(ik4), PARAMETER :: cb_XYL      = 19
 INTEGER(ik4), PARAMETER :: cb_ISOP     = 20
 INTEGER(ik4), PARAMETER :: cb_CRES     = 21
 INTEGER(ik4), PARAMETER :: cb_OPEN     = 22
 INTEGER(ik4), PARAMETER :: cb_ISOPN    = 23
 INTEGER(ik4), PARAMETER :: cb_ISOPRD   = 24
 INTEGER(ik4), PARAMETER :: cb_ONIT     = 25
 INTEGER(ik4), PARAMETER :: cb_MGLY     = 26
 INTEGER(ik4), PARAMETER :: cb_AONE     = 27
 INTEGER(ik4), PARAMETER :: cb_PAN      = 28
 INTEGER(ik4), PARAMETER :: cb_CH3OOH   = 29
 INTEGER(ik4), PARAMETER :: cb_ETHOOH   = 30
 INTEGER(ik4), PARAMETER :: cb_ALD2     = 31
 INTEGER(ik4), PARAMETER :: cb_HCHO     = 32
 INTEGER(ik4), PARAMETER :: cb_CH3OH    = 33

   real(rk8) , parameter :: w_no2 = 46.0D0
  real(rk8) , parameter :: w_no  = 30.0D0

  real(rk8) , parameter :: w_hono = 47.0D0
  real(rk8) , parameter :: w_hno2 = 47.0D0
  real(rk8) , parameter :: w_no3  = 62.0D0
  real(rk8) , parameter :: w_n2o5 = 108.0D0
  real(rk8) , parameter :: w_hno4 = 79.0D0
  real(rk8) , parameter :: w_hno3 = 63.0D0
  real(rk8) , parameter :: w_o3   = 48.0D0
  real(rk8) , parameter :: w_h2o2 = 34.0D0
  real(rk8) , parameter :: w_h2o  = 18.0D0

  real(rk8) , parameter :: w_so2  = 64.0D0
  real(rk8) , parameter :: w_sulf = 98.0D0
  real(rk8) , parameter :: w_h2so4= 98.0D0
  real(rk8) , parameter :: w_co   = 28.0D0
  real(rk8) , parameter :: w_co2  = 44.0D0
  real(rk8) , parameter :: w_h2   = 2.0D0

  real(rk8) , parameter :: w_oh   = 17.0D0
  real(rk8) , parameter :: w_ho2  = 33.0D0
  real(rk8) , parameter :: w_ro2  = 47.0D0
  real(rk8) , parameter :: w_xo2  = 47.0D0
  real(rk8) , parameter :: w_ch4    = 16.0D0
  real(rk8) , parameter :: w_ethan  = 30.0D0
  real(rk8) , parameter :: w_c2h6   = 30.070D0
  ! assumed molecular wieght for PAR (CBMZ mechanism)
  real(rk8) , parameter :: w_par    = 44.0D0
  real(rk8) , parameter :: w_hc3    = 44.0D0
  real(rk8) , parameter :: w_c3h8   = 44.10D0
  real(rk8) , parameter :: w_hc5    = 72.0D0
  real(rk8) , parameter :: w_hc8    = 114.0D0
  real(rk8) , parameter :: w_alk4   = 58.120D0
  real(rk8) , parameter :: w_alk7   = 100.200D0

  ! Alkene species in RADM2 and CBMZ
  real(rk8) , parameter :: w_ethene = 28.0D0
  real(rk8) , parameter :: w_eth = 28.0D0
  real(rk8) , parameter :: w_ol2    = 28.0D0
  real(rk8) , parameter :: w_olt    = 42.0D0
  real(rk8) , parameter :: w_oli    = 56.0D0
  real(rk8) , parameter :: w_olet   = 42.0D0
  real(rk8) , parameter :: w_olei   = 56.0D0
  real(rk8) , parameter :: w_prpe   = 42.0D0
  real(rk8) , parameter :: w_bute   = 56.0D0
  real(rk8) , parameter :: w_isop   = 68.0D0

  ! Aromatic
  real(rk8) , parameter :: w_tolu  = 92.0D0
  real(rk8) , parameter :: w_tol   = 92.0D0
  real(rk8) , parameter :: w_csl   = 108.0D0
  real(rk8) , parameter :: w_cres  = 108.0D0
  real(rk8) , parameter :: w_xyle  = 106.0D0
  real(rk8) , parameter :: w_xyl   = 106.0D0
  real(rk8) , parameter :: w_benz    = 78.110D0

  ! Carbonyls
  real(rk8) , parameter :: w_hcho    = 30.0D0
  real(rk8) , parameter :: w_ald2    = 44.0D0
  real(rk8) , parameter :: w_ket     = 72.0D0
  real(rk8) , parameter :: w_aone    = 72.0D0
  real(rk8) , parameter :: w_gly     = 58.0D0
  real(rk8) , parameter :: w_mgly    = 72.0D0
  ! Organic Nitrate
  real(rk8) , parameter :: w_pan     = 121.0D0
  real(rk8) , parameter :: w_tpan    = 147.0D0
  real(rk8) , parameter :: w_onit    = 119.0D0

  ! Organic Acids
  real(rk8), parameter  :: w_hcooh     = 46.0D0   !Formic acid
  real(rk8), parameter  :: w_ch3cooh   = 60.0D0   !Acetic acid
  real(rk8),  parameter :: w_rcooh   = 59.1D0

  ! Alcohol
  real(rk8), parameter  :: w_moh         = 32.0D0 !Methanol
  real(rk8), parameter  :: w_ch3oh       = 32.0D0 !Methanol
  real(rk8), parameter  :: w_eoh         = 46.0D0 !Ethanol
  real(rk8), parameter  :: w_c2h5oh      = 46.0D0 !Ethanol

  ! Organic Peroxid
  real(rk8), parameter  :: w_ch3ooh   = 48.0D0
  real(rk8), parameter  :: w_ethooh   = 74.0D0
  real(rk8) , parameter :: w_rooh    = 48.0D0
  ! Other species
  real(rk8) , parameter :: w_dms     = 62.0D0
  real(rk8) , parameter :: w_msa     = 96.0D0
  real(rk8) , parameter :: w_nh3     = 17.0D0
  real(rk8) , parameter :: w_apin    = 136.230D0
  real(rk8) , parameter :: w_limo    = 136.230D0

  ! intermediate species that do not undergo other process than chemistryi
  ! are assigned an arbitrary molecular weight , as anyway chemistry works
  ! with mol. and the conversion to mass is done for compatibility with
  ! other processes than chem.
  ! if these species are outputed in mass, the unit will be however wrong !

  real(rk8) , parameter :: w_O1D = 16.D0
  real(rk8) , parameter :: w_cro =48.D0
  real(rk8) , parameter :: w_to2 = 32.D0
  real(rk8) , parameter :: w_dummy = 1.D0
  real(rk8) , parameter :: w_open = 1.D0
  real(rk8) , parameter :: w_O3P = 48.D0
  real(rk8) , parameter :: w_isopn   = 68.0D0
  real(rk8) , parameter :: w_isopp   = 68.0D0
  real(rk8) , parameter :: w_isopo2   = 68.0D0
  real(rk8) , parameter :: w_isoprd   = 68.0D0
  real(rk8) , parameter :: w_ethp   = 28.0D0
  real(rk8) , parameter :: w_nap   = 1.0D0
  real(rk8) , parameter :: w_ch3o2   = 47.0D0
  real(rk8) , parameter :: w_ano2   =  46.0D0
  real(rk8) , parameter :: w_c2o3   =  72.0D0


end module mod_ch_param
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
