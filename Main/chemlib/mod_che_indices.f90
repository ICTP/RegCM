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

module mod_che_indices

  public

  !ah   gas_phase chemistry indecies for potential transported species 
  !ah   some of them correspond to the same compounds like isulf, ih2so4
  !ah   we did this just for flexibility
  integer :: iso2 , iso4 , idms , ibchl , ibchb , iochl , iochb

  integer :: imsa
  integer :: io3 , ino , ino2 , ino3 , ioh , iho2 , ih2o2
  integer :: ihno2 , ihno3 , ihno4
  integer :: isulf , ih2so4 , ihono , in2o5 , ihc , ihcr , ic2h4
  integer :: ico , ihcho , iald2 , iethe , ic2h6 , ic3h8
  integer :: iisop , itolue , ixyl , inh3 , ipan , in2o
  integer :: irooh , iacet , ibenz , ich4 , ico2
  integer :: inox , ihox , isox , ieoh , imoh , iaco2 , ircooh
  integer :: ipar , iolt , ioli
  integer :: iapin , ilimo
  integer :: ialk4, ialk7

  !*** abt added from wetdep scheme
  integer :: iisopno3 , ich3ooh , ihydrald , ihyac , ipooh
  integer :: ic3h7ooh , ic2h5ooh 
  integer :: iisopooh , imacrooh , ipb , ionit , ich3coooh , ich3cocho , ixooh
  integer :: ionitr , iglyald , imvk , imacr , isoa , inh4 , inh4no3 , ich3cooh
  integer :: iterpooh , itolooh , imekooh , ialkooh

  ! This the compound indecies
  ! this list for cbmz-sill mechanism
  integer , parameter :: ind_o3     = 1
  integer , parameter :: ind_no2    = 2 
  integer , parameter :: ind_no     = 3 
  integer , parameter :: ind_no3    = 4 
  integer , parameter :: ind_n2o5   = 5 
  integer , parameter :: ind_hno3   = 6 
  integer , parameter :: ind_hono   = 7 
  integer , parameter :: ind_hno4   = 8 
  integer , parameter :: ind_oh     = 9 
  integer , parameter :: ind_ho2    = 10 
  integer , parameter :: ind_h2o2   = 11 
  integer , parameter :: ind_so2    = 12 
  integer , parameter :: ind_sulf   = 13 
  integer , parameter :: ind_h2o    = 14 
  integer , parameter :: ind_co     = 15 
  integer , parameter :: ind_h2     = 16 
  integer , parameter :: ind_ch4    = 17 
  integer , parameter :: ind_c2h6   = 18 
  integer , parameter :: ind_par    = 19 !=3alk3+4.5alk4+alk7
  integer , parameter :: ind_moh    = 20 
  integer , parameter :: ind_roh    = 21
  integer , parameter :: ind_hcho   = 22 
  integer , parameter :: ind_ald2   = 23 
  integer , parameter :: ind_acet   = 24 
  integer , parameter :: ind_aco2   = 25 
  integer , parameter :: ind_rcooh  = 26
  integer , parameter :: ind_ethe   = 27 
  integer , parameter :: ind_prpe   = 28 
  integer , parameter :: ind_bute   = 29 
  integer , parameter :: ind_isop   = 30 
  integer , parameter :: ind_iprd   = 31 
  integer , parameter :: ind_tolu   = 32 
  integer , parameter :: ind_xyle   = 33 
  integer , parameter :: ind_ch3ooh = 34 
  integer , parameter :: ind_ethooh = 35 
  integer , parameter :: ind_rooh   = 36 
  integer , parameter :: ind_mgly   = 37
  integer , parameter :: ind_pan    = 38 
  integer , parameter :: ind_onit   = 39 
  integer , parameter :: ind_cres   = 40
  integer , parameter :: ind_cro    = 41
  integer , parameter :: ind_dms    = 42
  integer , parameter :: ind_open   = 43
  integer , parameter :: ind_prod   = 44
  integer , parameter :: ind_ch3o2  = 45
  integer , parameter :: ind_ethp   = 46
  integer , parameter :: ind_c2o3   = 47
  integer , parameter :: ind_ro2    = 48
  integer , parameter :: ind_to2    = 49
  integer , parameter :: ind_ano2   = 50
  integer , parameter :: ind_nap    = 51
  integer , parameter :: ind_isopp  = 52
  integer , parameter :: ind_isopn  = 53
  integer , parameter :: ind_isopo2 = 54
  integer , parameter :: ind_xo2    = 55
  integer , parameter :: ind_xpar   = 56
  integer , parameter :: ind_co2    = 57
  integer , parameter :: ind_nh3    = 58
  integer , parameter :: ind_hclg   = 59
  integer , parameter :: ind_naog   = 60
  integer , parameter :: ind_apin   = 61
  integer , parameter :: ind_limo   = 62
  integer , parameter :: ind_pip    = 63
  integer , parameter :: ind_pint   = 64
  integer , parameter :: ind_pio2   = 65
  integer , parameter :: ind_lio2   = 66
  integer , parameter :: ind_lip    = 67
  integer , parameter :: ind_pin2   = 68

end module mod_che_indices
