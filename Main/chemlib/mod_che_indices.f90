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

  use mod_intkinds

  public

  !ah   gas_phase chemistry indecies for potential transported species 
  !ah   some of them correspond to the same compounds like isulf, ih2so4
  !ah   we did this just for flexibility
  integer(ik4) :: iso2 , iso4 , idms , ibchl , ibchb , iochl , iochb

  integer(ik4) :: imsa
  integer(ik4) :: io3 , ino , ino2 , ino3 , ioh , iho2 , ih2o2
  integer(ik4) :: ihno2 , ihno3 , ihno4
  integer(ik4) :: isulf , ih2so4 , ihono , in2o5 , ihc , ihcr , ic2h4
  integer(ik4) :: ico , ihcho , iald2 , iethe , ic2h6 , ic3h8
  integer(ik4) :: iisop , itolue , ixyl , inh3 , ipan , in2o
  integer(ik4) :: irooh , iacet , ibenz , ich4 , ico2
  integer(ik4) :: inox , ihox , isox , ieoh , imoh , iaco2 , ircooh,ihcooh
  integer(ik4) :: ipar , iolt , ioli, imgly,icres,iopen,iisoprd, iethooh, ixo2
  integer(ik4) :: iapin , ilimo
  integer(ik4) :: ialk4, ialk7

  !*** abt added from wetdep scheme
  integer(ik4) :: iisopno3 , ich3ooh , ihydrald , ihyac , ipooh
  integer(ik4) :: ic3h7ooh , ic2h5ooh 
  integer(ik4) :: iisopooh , imacrooh , ipb , ionit , ich3coooh , ich3cocho , ixooh
  integer(ik4) :: ionitr , iglyald , imvk , imacr , isoa , inh4 , inh4no3 , ich3cooh
  integer(ik4) :: iterpooh , itolooh , imekooh , ialkooh

  ! This the compound indecies
  ! this list ishould absolutly be consistent with REACTION.DAT for Sillman/CBMZ mechanism
  integer(ik4) , parameter :: ind_o3     = 1
  integer(ik4) , parameter :: ind_no2    = 2 
  integer(ik4) , parameter :: ind_no     = 3 
  integer(ik4) , parameter :: ind_no3    = 4 
  integer(ik4) , parameter :: ind_n2o5   = 5 
  integer(ik4) , parameter :: ind_hno3   = 6 
  integer(ik4) , parameter :: ind_hono   = 7 
  integer(ik4) , parameter :: ind_hno4   = 8 
  integer(ik4) , parameter :: ind_oh     = 9 
  integer(ik4) , parameter :: ind_ho2    = 10 
  integer(ik4) , parameter :: ind_h2o2   = 11 
  integer(ik4) , parameter :: ind_so2    = 12 
  integer(ik4) , parameter :: ind_sulf   = 13 
  integer(ik4) , parameter :: ind_h2o    = 14 
  integer(ik4) , parameter :: ind_co     = 15 
  integer(ik4) , parameter :: ind_h2     = 16 
  integer(ik4) , parameter :: ind_ch4    = 17 
  integer(ik4) , parameter :: ind_c2h6   = 18 
  integer(ik4) , parameter :: ind_par    = 19 !=3alk3+4.5alk4+alk7
  integer(ik4) , parameter :: ind_moh    = 20 
  integer(ik4) , parameter :: ind_roh    = 21
  integer(ik4) , parameter :: ind_hcho   = 22 
  integer(ik4) , parameter :: ind_ald2   = 23 
  integer(ik4) , parameter :: ind_acet   = 24 
  integer(ik4) , parameter :: ind_hcooh   = 25 
  integer(ik4) , parameter :: ind_rcooh  = 26
  integer(ik4) , parameter :: ind_ethe   = 27 
  integer(ik4) , parameter :: ind_prpe   = 28 
  integer(ik4) , parameter :: ind_bute   = 29 
  integer(ik4) , parameter :: ind_isop   = 30 
  integer(ik4) , parameter :: ind_isoprd   = 31 
  integer(ik4) , parameter :: ind_tolu   = 32 
  integer(ik4) , parameter :: ind_xyle   = 33 
  integer(ik4) , parameter :: ind_ch3ooh = 34 
  integer(ik4) , parameter :: ind_ethooh = 35 
  integer(ik4) , parameter :: ind_rooh   = 36 
  integer(ik4) , parameter :: ind_mgly   = 37
  integer(ik4) , parameter :: ind_pan    = 38 
  integer(ik4) , parameter :: ind_onit   = 39 
  integer(ik4) , parameter :: ind_cres   = 40
  integer(ik4) , parameter :: ind_cro    = 41
  integer(ik4) , parameter :: ind_dms    = 42
  integer(ik4) , parameter :: ind_open   = 43
  integer(ik4) , parameter :: ind_prod   = 44
  integer(ik4) , parameter :: ind_ch3o2  = 45
  integer(ik4) , parameter :: ind_ethp   = 46
  integer(ik4) , parameter :: ind_c2o3   = 47
  integer(ik4) , parameter :: ind_ro2    = 48
  integer(ik4) , parameter :: ind_to2    = 49
  integer(ik4) , parameter :: ind_ano2   = 50
  integer(ik4) , parameter :: ind_nap    = 51
  integer(ik4) , parameter :: ind_isopp  = 52
  integer(ik4) , parameter :: ind_isopn  = 53
  integer(ik4) , parameter :: ind_isopo2 = 54
  integer(ik4) , parameter :: ind_xo2    = 55
  integer(ik4) , parameter :: ind_xpar   = 56
  integer(ik4) , parameter :: ind_co2    = 57
  integer(ik4) , parameter :: ind_nh3    = 58
  integer(ik4) , parameter :: ind_hclg   = 59
  integer(ik4) , parameter :: ind_naog   = 60
  integer(ik4) , parameter :: ind_apin   = 61
  integer(ik4) , parameter :: ind_limo   = 62
  integer(ik4) , parameter :: ind_pip    = 63
  integer(ik4) , parameter :: ind_pint   = 64
  integer(ik4) , parameter :: ind_pio2   = 65
  integer(ik4) , parameter :: ind_lio2   = 66
  integer(ik4) , parameter :: ind_lip    = 67
  integer(ik4) , parameter :: ind_pin2   = 68

! indices for jvalues 
  integer(ik4) , parameter :: jvO2 = 1
  integer(ik4) , parameter :: jvO3a = 2
  integer(ik4) , parameter :: jvO3b = 3
  integer(ik4) , parameter :: jvNO2 = 4
  integer(ik4) , parameter :: jvNO3a = 5
  integer(ik4) , parameter :: jvNO3b = 6
  integer(ik4) , parameter :: jvN2O5a = 7
  integer(ik4) , parameter :: jvN2O5b = 8
  integer(ik4) , parameter :: jvN2O = 9
  integer(ik4) , parameter :: jvHO2 = 10
  integer(ik4) , parameter :: jvH2O2 = 11
  integer(ik4) , parameter :: jvHNO2 = 12
  integer(ik4) , parameter :: jvHNO3 = 13
  integer(ik4) , parameter :: jvHNO4 = 14
  integer(ik4) , parameter :: jvCH2Oa = 15
  integer(ik4) , parameter :: jvCH2Ob = 16
  integer(ik4) , parameter :: jvCH3CHOa = 17
  integer(ik4) , parameter :: jvCH3CHOb = 18
  integer(ik4) , parameter :: jvCH3CHOc = 19
  integer(ik4) , parameter :: jvC2H5CHO = 20
  integer(ik4) , parameter :: jvCHOCHO = 21
  integer(ik4) , parameter :: jvCH3COCHO = 22
  integer(ik4) , parameter :: jvCH3COCH3 = 23
  integer(ik4) , parameter :: jvCH3OOH = 24
  integer(ik4) , parameter :: jvCH3ONO2 = 25
  integer(ik4) , parameter :: jvPAN = 26

end module mod_che_indices
