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
  integer(ik4) , parameter :: mz_NO      = 1
  integer(ik4) , parameter :: mz_NO2     = 2
  integer(ik4) , parameter :: mz_N2O5    = 3
  integer(ik4) , parameter :: mz_HNO3    = 4
  integer(ik4) , parameter :: mz_HO2NO2  = 5
  integer(ik4) , parameter :: mz_O3      = 6
  integer(ik4) , parameter :: mz_H2O2    = 7
  integer(ik4) , parameter :: mz_SO2     = 8
  integer(ik4) , parameter :: mz_SO4     = 9
  integer(ik4) , parameter :: mz_CH4     = 10
  integer(ik4) , parameter :: mz_CH2O    = 11
  integer(ik4) , parameter :: mz_CH3OH   = 12
  integer(ik4) , parameter :: mz_PAN     = 13
  integer(ik4) , parameter :: mz_C2H6    = 14
  integer(ik4) , parameter :: mz_C3H8    = 15
  integer(ik4) , parameter :: mz_BIGALK  = 16
  integer(ik4) , parameter :: mz_C2H4    = 17
  integer(ik4) , parameter :: mz_C3H6    = 18
  integer(ik4) , parameter :: mz_BIGENE  = 19
  integer(ik4) , parameter :: mz_TOLUENE = 20
  integer(ik4) , parameter :: mz_ISOP    = 21
  integer(ik4) , parameter :: mz_CH3CHO  = 22
  integer(ik4) , parameter :: mz_CH3COOH = 23
  integer(ik4) , parameter :: mz_GLYALD  = 24
  integer(ik4) , parameter :: mz_CH3OOH  = 25
  integer(ik4) , parameter :: mz_C2H5OOH = 26
  integer(ik4) , parameter :: mz_CH3COCH3= 27
  integer(ik4) , parameter :: mz_HYAC    = 28
  integer(ik4) , parameter :: mz_CH3COCHO =29
  integer(ik4) , parameter :: mz_ONIT    = 30
  integer(ik4) , parameter :: mz_MEK     = 31
  integer(ik4) , parameter :: mz_MVK     = 32
  integer(ik4) , parameter :: mz_MACR    = 33
  integer(ik4) , parameter :: mz_HYDRALD = 34
  integer(ik4) , parameter :: mz_BIGALD  = 35
  integer(ik4) , parameter :: mz_ISOPNO3 = 36
  integer(ik4) , parameter :: mz_ONITR   = 37
  integer(ik4) , parameter :: mz_CRESOL  = 38
  integer(ik4) , parameter :: mz_CO      = 39
  integer(ik4) , parameter :: mz_DMS     = 40

  !CBMZ indcies
  integer(ik4) , parameter :: cb_O3       = 1
  integer(ik4) , parameter :: cb_NO       = 2
  integer(ik4) , parameter :: cb_NO2      = 3
  integer(ik4) , parameter :: cb_HNO3     = 4
  integer(ik4) , parameter :: cb_HNO4     = 5
  integer(ik4) , parameter :: cb_N2O5     = 6
  integer(ik4) , parameter :: cb_H2O2     = 7
  integer(ik4) , parameter :: cb_CH4      = 8
  integer(ik4) , parameter :: cb_CO       = 9
  integer(ik4) , parameter :: cb_SO2      = 10
  integer(ik4) , parameter :: cb_H2SO4    = 11
  integer(ik4) , parameter :: cb_DMS      = 12
  integer(ik4) , parameter :: cb_PAR      = 13
  integer(ik4) , parameter :: cb_C2H6     = 14
  integer(ik4) , parameter :: cb_ETH      = 15
  integer(ik4) , parameter :: cb_OLET     = 16
  integer(ik4) , parameter :: cb_OLEI     = 17
  integer(ik4) , parameter :: cb_TOL      = 18
  integer(ik4) , parameter :: cb_XYL      = 19
  integer(ik4) , parameter :: cb_ISOP     = 20
  integer(ik4) , parameter :: cb_CRES     = 21
  integer(ik4) , parameter :: cb_OPEN     = 22
  integer(ik4) , parameter :: cb_ISOPN    = 23
  integer(ik4) , parameter :: cb_ISOPRD   = 24
  integer(ik4) , parameter :: cb_ONIT     = 25
  integer(ik4) , parameter :: cb_MGLY     = 26
  integer(ik4) , parameter :: cb_AONE     = 27
  integer(ik4) , parameter :: cb_PAN      = 28
  integer(ik4) , parameter :: cb_CH3OOH   = 29
  integer(ik4) , parameter :: cb_ETHOOH   = 30
  integer(ik4) , parameter :: cb_ALD2     = 31
  integer(ik4) , parameter :: cb_HCHO     = 32
  integer(ik4) , parameter :: cb_CH3OH    = 33

  real(rkx) , parameter :: w_no2 = 46.0_rkx
  real(rkx) , parameter :: w_no  = 30.0_rkx

  real(rkx) , parameter :: w_hono = 47.0_rkx
  real(rkx) , parameter :: w_hno2 = 47.0_rkx
  real(rkx) , parameter :: w_no3  = 62.0_rkx
  real(rkx) , parameter :: w_n2o5 = 108.0_rkx
  real(rkx) , parameter :: w_hno4 = 79.0_rkx
  real(rkx) , parameter :: w_hno3 = 63.0_rkx
  real(rkx) , parameter :: w_o3   = 48.0_rkx
  real(rkx) , parameter :: w_h2o2 = 34.0_rkx
  real(rkx) , parameter :: w_h2o  = 18.0_rkx

  real(rkx) , parameter :: w_so2  = 64.0_rkx
  real(rkx) , parameter :: w_sulf = 98.0_rkx
  real(rkx) , parameter :: w_h2so4= 98.0_rkx
  real(rkx) , parameter :: w_co   = 28.0_rkx
  real(rkx) , parameter :: w_co2  = 44.0_rkx
  real(rkx) , parameter :: w_h2   = 2.0_rkx

  real(rkx) , parameter :: w_oh   = 17.0_rkx
  real(rkx) , parameter :: w_ho2  = 33.0_rkx
  real(rkx) , parameter :: w_ro2  = 47.0_rkx
  real(rkx) , parameter :: w_xo2  = 47.0_rkx
  real(rkx) , parameter :: w_ch4    = 16.0_rkx
  real(rkx) , parameter :: w_ethan  = 30.0_rkx
  real(rkx) , parameter :: w_c2h6   = 30.070_rkx
  ! assumed molecular wieght for PAR (CBMZ mechanism)
  real(rkx) , parameter :: w_par    = 44.0_rkx
  real(rkx) , parameter :: w_hc3    = 44.0_rkx
  real(rkx) , parameter :: w_c3h8   = 44.10_rkx
  real(rkx) , parameter :: w_hc5    = 72.0_rkx
  real(rkx) , parameter :: w_hc8    = 114.0_rkx
  real(rkx) , parameter :: w_alk4   = 58.120_rkx
  real(rkx) , parameter :: w_alk7   = 100.200_rkx

  ! Alkene species in RADM2 and CBMZ
  real(rkx) , parameter :: w_ethene = 28.0_rkx
  real(rkx) , parameter :: w_eth = 28.0_rkx
  real(rkx) , parameter :: w_ol2    = 28.0_rkx
  real(rkx) , parameter :: w_olt    = 42.0_rkx
  real(rkx) , parameter :: w_oli    = 56.0_rkx
  real(rkx) , parameter :: w_olet   = 42.0_rkx
  real(rkx) , parameter :: w_olei   = 56.0_rkx
  real(rkx) , parameter :: w_prpe   = 42.0_rkx
  real(rkx) , parameter :: w_bute   = 56.0_rkx
  real(rkx) , parameter :: w_isop   = 68.0_rkx

  ! Aromatic
  real(rkx) , parameter :: w_tolu  = 92.0_rkx
  real(rkx) , parameter :: w_tol   = 92.0_rkx
  real(rkx) , parameter :: w_csl   = 108.0_rkx
  real(rkx) , parameter :: w_cres  = 108.0_rkx
  real(rkx) , parameter :: w_xyle  = 106.0_rkx
  real(rkx) , parameter :: w_xyl   = 106.0_rkx
  real(rkx) , parameter :: w_benz    = 78.110_rkx

  ! Carbonyls
  real(rkx) , parameter :: w_hcho    = 30.0_rkx
  real(rkx) , parameter :: w_ald2    = 44.0_rkx
  real(rkx) , parameter :: w_ket     = 72.0_rkx
  real(rkx) , parameter :: w_aone    = 72.0_rkx
  real(rkx) , parameter :: w_gly     = 58.0_rkx
  real(rkx) , parameter :: w_mgly    = 72.0_rkx
  ! Organic Nitrate
  real(rkx) , parameter :: w_pan     = 121.0_rkx
  real(rkx) , parameter :: w_tpan    = 147.0_rkx
  real(rkx) , parameter :: w_onit    = 119.0_rkx

  ! Organic Acids
  real(rkx), parameter  :: w_hcooh     = 46.0_rkx   !Formic acid
  real(rkx), parameter  :: w_ch3cooh   = 60.0_rkx   !Acetic acid
  real(rkx),  parameter :: w_rcooh   = 59.1_rkx

  ! Alcohol
  real(rkx), parameter  :: w_moh         = 32.0_rkx !Methanol
  real(rkx), parameter  :: w_ch3oh       = 32.0_rkx !Methanol
  real(rkx), parameter  :: w_eoh         = 46.0_rkx !Ethanol
  real(rkx), parameter  :: w_c2h5oh      = 46.0_rkx !Ethanol

  ! Organic Peroxid
  real(rkx), parameter  :: w_ch3ooh   = 48.0_rkx
  real(rkx), parameter  :: w_ethooh   = 74.0_rkx
  real(rkx) , parameter :: w_rooh    = 48.0_rkx
  ! Other species
  real(rkx) , parameter :: w_dms     = 62.0_rkx
  real(rkx) , parameter :: w_msa     = 96.0_rkx
  real(rkx) , parameter :: w_nh3     = 17.0_rkx
  real(rkx) , parameter :: w_apin    = 136.230_rkx
  real(rkx) , parameter :: w_limo    = 136.230_rkx

  ! intermediate species that do not undergo other process than chemistryi
  ! are assigned an arbitrary molecular weight , as anyway chemistry works
  ! with mol. and the conversion to mass is done for compatibility with
  ! other processes than chem.
  ! if these species are outputed in mass, the unit will be however wrong !

  real(rkx) , parameter :: w_O1D = 16._rkx
  real(rkx) , parameter :: w_cro =48._rkx
  real(rkx) , parameter :: w_to2 = 32._rkx
  real(rkx) , parameter :: w_dummy = 1._rkx
  real(rkx) , parameter :: w_open = 1._rkx
  real(rkx) , parameter :: w_O3P = 48._rkx
  real(rkx) , parameter :: w_isopn   = 68.0_rkx
  real(rkx) , parameter :: w_isopp   = 68.0_rkx
  real(rkx) , parameter :: w_isopo2   = 68.0_rkx
  real(rkx) , parameter :: w_isoprd   = 68.0_rkx
  real(rkx) , parameter :: w_ethp   = 28.0_rkx
  real(rkx) , parameter :: w_nap   = 1.0_rkx
  real(rkx) , parameter :: w_ch3o2   = 47.0_rkx
  real(rkx) , parameter :: w_ano2   =  46.0_rkx
  real(rkx) , parameter :: w_c2o3   =  72.0_rkx

end module mod_ch_param
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
