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
  integer(ik4) :: ipar , iolt , ioli, imgly,icres,iopen,iisoprd, iethooh, ixo2,iro2
  integer(ik4) :: iapin , ilimo
  integer(ik4) :: ialk4, ialk7
  integer(ik4) :: ianh4, iano3
  

  !*** abt added from wetdep scheme
  integer(ik4) :: iisopno3 , ich3ooh , ihydrald , ihyac , ipooh
  integer(ik4) :: ic3h7ooh , ic2h5ooh 
  integer(ik4) :: iisopooh , imacrooh , ipb , ionit , ich3coooh , ich3cocho , ixooh
  integer(ik4) :: ionitr , iglyald , imvk , imacr , isoa , inh4 , inh4no3 , ich3cooh
  integer(ik4) :: iterpooh , itolooh , imekooh , ialkooh

  !
  integer(ik4) :: ipollen

  ! This the compound indecies
  ! this list ishould absolutly be consistent with KPP GAS_CBMZ/mod_cbmz_Parameters.f90 file
  INTEGER(ik4), PARAMETER :: ind_CO2 = 1
  INTEGER(ik4), PARAMETER :: ind_H2SO4 = 2
  INTEGER(ik4), PARAMETER :: ind_HCOOH = 3
  INTEGER(ik4), PARAMETER :: ind_RCOOH = 4
  INTEGER(ik4), PARAMETER :: ind_MSA = 5
  INTEGER(ik4), PARAMETER :: ind_DUMMY = 6
  INTEGER(ik4), PARAMETER :: ind_PAN = 7
  INTEGER(ik4), PARAMETER :: ind_TOL = 8
  INTEGER(ik4), PARAMETER :: ind_O1D = 9
  INTEGER(ik4), PARAMETER :: ind_H2O2 = 10
  INTEGER(ik4), PARAMETER :: ind_SO2 = 11
  INTEGER(ik4), PARAMETER :: ind_XYL = 12
  INTEGER(ik4), PARAMETER :: ind_CH4 = 13
  INTEGER(ik4), PARAMETER :: ind_C2H6 = 14
  INTEGER(ik4), PARAMETER :: ind_CRO = 15
  INTEGER(ik4), PARAMETER :: ind_DMS = 16
  INTEGER(ik4), PARAMETER :: ind_HNO4 = 17
  INTEGER(ik4), PARAMETER :: ind_H2 = 18
  INTEGER(ik4), PARAMETER :: ind_TO2 = 19
  INTEGER(ik4), PARAMETER :: ind_CH3OH = 20
  INTEGER(ik4), PARAMETER :: ind_HNO2 = 21
  INTEGER(ik4), PARAMETER :: ind_CH3OOH = 22
  INTEGER(ik4), PARAMETER :: ind_ETHOOH = 23
  INTEGER(ik4), PARAMETER :: ind_N2O5 = 24
  INTEGER(ik4), PARAMETER :: ind_ETH = 25
  INTEGER(ik4), PARAMETER :: ind_CRES = 26
  INTEGER(ik4), PARAMETER :: ind_O3P = 27
  INTEGER(ik4), PARAMETER :: ind_CO = 28
  INTEGER(ik4), PARAMETER :: ind_HNO3 = 29
  INTEGER(ik4), PARAMETER :: ind_PAR = 30
  INTEGER(ik4), PARAMETER :: ind_OPEN = 31
  INTEGER(ik4), PARAMETER :: ind_ISOPN = 32
  INTEGER(ik4), PARAMETER :: ind_ISOPP = 33
  INTEGER(ik4), PARAMETER :: ind_ISOPO2 = 34
  INTEGER(ik4), PARAMETER :: ind_H2O = 35
  INTEGER(ik4), PARAMETER :: ind_AONE = 36
  INTEGER(ik4), PARAMETER :: ind_OLEI = 37
  INTEGER(ik4), PARAMETER :: ind_ISOP = 38
  INTEGER(ik4), PARAMETER :: ind_HCHO = 39
  INTEGER(ik4), PARAMETER :: ind_OLET = 40
  INTEGER(ik4), PARAMETER :: ind_XO2 = 41
  INTEGER(ik4), PARAMETER :: ind_MGLY = 42
  INTEGER(ik4), PARAMETER :: ind_ETHP = 43
  INTEGER(ik4), PARAMETER :: ind_NAP = 44
  INTEGER(ik4), PARAMETER :: ind_ALD2 = 45
  INTEGER(ik4), PARAMETER :: ind_CH3O2 = 46
  INTEGER(ik4), PARAMETER :: ind_ISOPRD = 47
  INTEGER(ik4), PARAMETER :: ind_ANO2 = 48
  INTEGER(ik4), PARAMETER :: ind_ROOH = 49
  INTEGER(ik4), PARAMETER :: ind_RO2 = 50
  INTEGER(ik4), PARAMETER :: ind_ONIT = 51
  INTEGER(ik4), PARAMETER :: ind_HO2 = 52
  INTEGER(ik4), PARAMETER :: ind_O3 = 53
  INTEGER(ik4), PARAMETER :: ind_OH = 54
  INTEGER(ik4), PARAMETER :: ind_NO = 55
  INTEGER(ik4), PARAMETER :: ind_NO2 = 56
  INTEGER(ik4), PARAMETER :: ind_NO3 = 57
  INTEGER(ik4), PARAMETER :: ind_C2O3 = 58

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
