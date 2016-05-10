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

module mod_che_molwg

  use mod_intkinds
  use mod_realkinds
  use mod_che_indices

  implicit none

  public

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

  ! Alkane species in RADM2 and CBMZ
  ! some species are repeate with differen name convention
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

  ! define here a table of molecular weight for the CBMZ species.

  real(rkx),dimension(totsp) ::  mw_cbmz

  data mw_cbmz / W_CO2 ,   & ! 1
                 W_H2SO4,  & ! 2
                 W_HCOOH,  & ! 3
                 W_RCOOH,  & ! 4
                 W_MSA,    & ! 5
                 W_DUMMY,  & ! 6
                 W_PAN,    & ! 7
                 W_TOL,    & ! 8
                 W_O1D,    & ! 9
                 W_H2O2,   & ! 10
                 W_SO2,    & ! 11
                 W_XYL,    & ! 12
                 W_CH4,    & ! 13
                 W_C2H6,   & ! 14
                 W_CRO,    & ! 15
                 W_DMS,    & ! 16
                 W_HNO4,   & ! 17
                 W_H2,     & ! 18
                 W_TO2,    & ! 19
                 W_CH3OH,  & ! 20
                 W_HNO2,   & ! 21
                 W_CH3OOH, & ! 22
                 W_ETHOOH, & ! 23
                 W_N2O5,   & ! 24
                 W_ETH,    & ! 25
                 W_CRES,   & ! 26
                 W_O3P,    & ! 27
                 W_CO,     & ! 28
                 W_HNO3,   & ! 29
                 W_PAR,    & ! 30
                 W_OPEN,   & ! 31
                 W_ISOPN,  & ! 32
                 W_ISOPP,  & ! 33
                 W_ISOPO2, & ! 34
                 W_H2O,    & ! 35
                 W_AONE,   & ! 36
                 W_OLEI,   & ! 37
                 W_ISOP,   & ! 38
                 W_HCHO,   & ! 39
                 W_OLET,   & ! 40
                 W_XO2,    & ! 41
                 W_MGLY,   & ! 42
                 W_ETHP,   & ! 43
                 W_NAP,    & ! 44
                 W_ALD2,   & ! 45
                 W_CH3O2,  & ! 46
                 W_ISOPRD, & ! 47
                 W_ANO2,   & ! 48
                 W_ROOH,   & ! 49
                 W_RO2,    & ! 50
                 W_ONIT,   & ! 51
                 W_HO2,    & ! 52
                 W_O3,     & ! 53
                 W_OH,     & ! 54
                 W_NO,     & ! 55
                 W_NO2,    & ! 56
                 W_NO3,    & ! 57
                 W_C2O3 /    ! 58

end module mod_che_molwg

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
