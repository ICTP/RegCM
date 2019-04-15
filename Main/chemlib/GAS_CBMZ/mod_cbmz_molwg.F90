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

module mod_cbmz_molwg
!
  use mod_intkinds
  use mod_realkinds
!
  public
!
  real(rkx) , parameter :: w_no2 = 46.0_rkx
  real(rkx) , parameter :: w_no  = 30.0_rkx

  real(rkx) , parameter :: w_hono = 47.0_rkx
  real(rkx) , parameter :: w_no3  = 62.0_rkx
  real(rkx) , parameter :: w_n2o5 = 108.0_rkx
  real(rkx) , parameter :: w_hno4 = 79.0_rkx
  real(rkx) , parameter :: w_hno3 = 63.0_rkx
  real(rkx) , parameter :: w_o3   = 48.0_rkx
  real(rkx) , parameter :: w_h2o2 = 34.0_rkx

  real(rkx) , parameter :: w_so2  = 64.0_rkx
  real(rkx) , parameter :: w_sulf = 98.0_rkx
  real(rkx) , parameter :: w_co   = 28.0_rkx
  real(rkx) , parameter :: w_co2  = 44.0_rkx
  real(rkx) , parameter :: w_h2   = 2.0_rkx

  real(rkx) , parameter :: w_oh   = 17.0_rkx
  real(rkx) , parameter :: w_ho2  = 33.0_rkx

  real(rkx) , parameter :: w_ch4    = 16.0_rkx
  real(rkx) , parameter :: w_ethan  = 30.0_rkx
  real(rkx) , parameter :: w_hc3    = 44.0_rkx
  real(rkx) , parameter :: w_hc5    = 72.0_rkx
  real(rkx) , parameter :: w_hc8    = 114.0_rkx

  real(rkx) , parameter :: w_ethene = 28.0_rkx
  real(rkx) , parameter :: w_ol2    = 28.0_rkx
  real(rkx) , parameter :: w_olt    = 42.0_rkx
  real(rkx) , parameter :: w_prpe   = 42.0_rkx
  real(rkx) , parameter :: w_bute   = 56.0_rkx
  real(rkx) , parameter :: w_oli    = 56.0_rkx
  real(rkx) , parameter :: w_isop   = 68.0_rkx

  real(rkx) , parameter :: w_tolu   = 92.0_rkx
  real(rkx) , parameter :: w_csl    = 108.0_rkx
  real(rkx) , parameter :: w_xyle   = 106.0_rkx

  real(rkx) , parameter :: w_hcho    = 30.0_rkx
  real(rkx) , parameter :: w_ald2    = 44.0_rkx
  real(rkx) , parameter :: w_ket     = 72.0_rkx
  real(rkx) , parameter :: w_gly     = 58.0_rkx
  real(rkx) , parameter :: w_mgly    = 72.0_rkx

  real(rkx) , parameter :: w_pan     = 121.0_rkx
  real(rkx) , parameter :: w_tpan    = 147.0_rkx
  real(rkx) , parameter :: w_alco    = 32.0_rkx

  real(rkx), parameter  :: w_cres   = 108.0_rkx
  real(rkx), parameter  :: w_hcooh   = 46.0_rkx !
  real(rkx), parameter  :: w_ch3ooh   = 48.0_rkx !
  real(rkx), parameter  :: w_ethooh   = 74.0_rkx !


  real(rkx) , parameter :: w_dms     = 62.0_rkx
  real(rkx) , parameter :: w_rooh    = 48.0_rkx
  real(rkx) , parameter :: w_nh3     = 17.0_rkx
  real(rkx) , parameter :: w_c2h6    = 30.070_rkx
  real(rkx) , parameter :: w_c3h8    = 44.10_rkx
  real(rkx) , parameter :: w_alk4    = 58.120_rkx
  real(rkx) , parameter :: w_alk7    = 100.200_rkx
  real(rkx) , parameter :: w_mo2     = 47.0_rkx
  real(rkx) , parameter :: w_acet    = 58.080_rkx
  real(rkx) , parameter :: w_moh     = 32.040_rkx
  real(rkx) , parameter :: w_eoh     = 46.070_rkx
  real(rkx) , parameter :: w_benz    = 78.110_rkx
  real(rkx),  parameter :: w_rcooh   = 59.1_rkx
  real(rkx) , parameter :: w_apin    = 136.230_rkx
  real(rkx) , parameter :: w_limo    = 136.230_rkx
!
!here are some intermediate species or operator : set to molecular weight = 1, or determine form Emmons et al GMDD
!It should not matter since evry thing is converted to mol before the chem solver , so proprtion between reactants are
!conserved , unless other processes than chemical reaction modify concentrations
!BUT unit is wong in the output
  real(rkx) , parameter :: w_xo2    = 149.0_rkx
  real(rkx), parameter  :: w_open   = 1.0_rkx !
  real(rkx), parameter  :: w_isoprd =  116._rkx !
  real(rkx), parameter  :: w_onit   = 119.0_rkx !


end module mod_cbmz_molwg
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
