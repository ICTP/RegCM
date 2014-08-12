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

  implicit none

  public

  !ah   gas_phase chemistry indecies for potential transported species 
  !ah   some of them correspond to the same compounds like isulf, ih2so4
  !ah   we did this just for flexibility
 
  integer(ik4) :: iso2 , iso4 , idms , ibchl , ibchb , iochl , iochb

  integer(ik4) :: imsa
  integer(ik4) :: io3 , ino , ino2 , ino3 , ioh , iho2 , ih2o2
  integer(ik4) :: ihno2 , ihno3 , ihno4
  integer(ik4) :: isulf , ih2so4 , ihono , in2o5 , ihc , ihcr , ic2h4
  integer(ik4) :: ico , ihcho , iald2 , ieth , ic2h6 , ic3h8,ic3h6
  integer(ik4) :: iisop , itol , ixyl , inh3 , ipan , in2o
  integer(ik4) :: irooh , iaone , ibenz , ich4 , ico2
  integer(ik4) :: inox , ihox , isox , ieoh , ich3oh , iaco2 , ircooh,ihcooh
  integer(ik4) :: ipar , iolet , iolei , imgly , icres , iopen , iisoprd,iisopn
  integer(ik4) :: iethooh , ixo2 , iro2
  integer(ik4) :: iapin , ilimo
  integer(ik4) :: ialk4, ialk7
  integer(ik4) :: ianh4, iano3

  !*** abt added from wetdep scheme
  integer(ik4) :: iisopno3 , ich3ooh , ihydrald , ihyac , ipooh
  integer(ik4) :: ic3h7ooh , ic2h5ooh 
  integer(ik4) :: iisopooh , imacrooh , ipb , ionit , ich3coooh
  integer(ik4) :: ich3cocho , ixooh
  integer(ik4) :: ionitr , iglyald , imvk , imacr , isoa , inh4
  integer(ik4) :: inh4no3 , ich3cooh
  integer(ik4) :: iterpooh , itolooh , imekooh , ialkooh

  !
  integer(ik4) :: ipollen

  ! list and name of cbmz species : must be absolutely consistant with 
  ! mod_cbmz_Parameters 

  integer(ik4) , parameter :: totsp = 58
  character(len=8),target, dimension(totsp) :: cbmzspec
  data  cbmzspec /'CO2',     & ! 1
                  'H2SO4',   & ! 2
                  'HCOOH',   & ! 3
                  'RCOOH',   & ! 4
                  'MSA',     & ! 5
                  'DUMMY',   & ! 6
                  'PAN',     & ! 7
                  'TOL',     & ! 8
                  'O1D',     & ! 9
                  'H2O2',    & ! 10
                  'SO2',     & ! 11
                  'XYL',     & ! 12
                  'CH4',     & ! 13
                  'C2H6',    & ! 14
                  'CRO',     & ! 15
                  'DMS',     & ! 16
                  'HNO4',    & ! 17
                  'H2',      & ! 18
                  'TO2',     & ! 19
                  'CH3OH',   & ! 20
                  'HNO2',    & ! 21
                  'CH3OOH',  & ! 22
                  'ETHOOH',  & ! 23
                  'N2O5',    & ! 24
                  'ETH',     & ! 25
                  'CRES',    & ! 26
                  'O3P',     & ! 27
                  'CO',      & ! 28
                  'HNO3',    & ! 29
                  'PAR',     & ! 30
                  'OPEN',    & ! 31
                  'ISOPN',   & ! 32
                  'ISOPP',   & ! 33
                  'ISOPO2',  & ! 34
                  'H2O',     & ! 35
                  'AONE',    & ! 36
                  'OLEI',    & ! 37
                  'ISOP',    & ! 38
                  'HCHO',    & ! 39
                  'OLET',    & ! 40
                  'XO2',     & ! 41
                  'MGLY',    & ! 42
                  'ETHP',    & ! 43
                  'NAP',     & ! 44
                  'ALD2',    & ! 45
                  'CH3O2',   & ! 46
                  'ISOPRD',  & ! 47
                  'ANO2',    & ! 48
                  'ROOH',    & ! 49
                  'RO2',     & ! 50
                  'ONIT',    & ! 51
                  'HO2',     & ! 52
                  'O3',      & ! 53
                  'OH',      & ! 54
                  'NO',      & ! 55
                  'NO2',     & ! 56
                  'NO3',     & ! 57
                  'C2O3'/      ! 58

  integer(ik4) , parameter :: ind_CO2 = 1
  integer(ik4) , parameter :: ind_H2SO4 = 2
  integer(ik4) , parameter :: ind_HCOOH = 3
  integer(ik4) , parameter :: ind_RCOOH = 4
  integer(ik4) , parameter :: ind_MSA = 5
  integer(ik4) , parameter :: ind_DUMMY = 6
  integer(ik4) , parameter :: ind_PAN = 7
  integer(ik4) , parameter :: ind_TOL = 8
  integer(ik4) , parameter :: ind_O1D = 9
  integer(ik4) , parameter :: ind_H2O2 = 10
  integer(ik4) , parameter :: ind_SO2 = 11
  integer(ik4) , parameter :: ind_XYL = 12
  integer(ik4) , parameter :: ind_CH4 = 13
  integer(ik4) , parameter :: ind_C2H6 = 14
  integer(ik4) , parameter :: ind_CRO = 15
  integer(ik4) , parameter :: ind_DMS = 16
  integer(ik4) , parameter :: ind_HNO4 = 17
  integer(ik4) , parameter :: ind_H2 = 18
  integer(ik4) , parameter :: ind_TO2 = 19
  integer(ik4) , parameter :: ind_CH3OH = 20
  integer(ik4) , parameter :: ind_HNO2 = 21
  integer(ik4) , parameter :: ind_CH3OOH = 22
  integer(ik4) , parameter :: ind_ETHOOH = 23
  integer(ik4) , parameter :: ind_N2O5 = 24
  integer(ik4) , parameter :: ind_ETH = 25
  integer(ik4) , parameter :: ind_CRES = 26
  integer(ik4) , parameter :: ind_O3P = 27
  integer(ik4) , parameter :: ind_CO = 28
  integer(ik4) , parameter :: ind_HNO3 = 29
  integer(ik4) , parameter :: ind_PAR = 30
  integer(ik4) , parameter :: ind_OPEN = 31
  integer(ik4) , parameter :: ind_ISOPN = 32
  integer(ik4) , parameter :: ind_ISOPP = 33
  integer(ik4) , parameter :: ind_ISOPO2 = 34
  integer(ik4) , parameter :: ind_H2O = 35
  integer(ik4) , parameter :: ind_AONE = 36
  integer(ik4) , parameter :: ind_OLEI = 37
  integer(ik4) , parameter :: ind_ISOP = 38
  integer(ik4) , parameter :: ind_HCHO = 39
  integer(ik4) , parameter :: ind_OLET = 40
  integer(ik4) , parameter :: ind_XO2 = 41
  integer(ik4) , parameter :: ind_MGLY = 42
  integer(ik4) , parameter :: ind_ETHP = 43
  integer(ik4) , parameter :: ind_NAP = 44
  integer(ik4) , parameter :: ind_ALD2 = 45
  integer(ik4) , parameter :: ind_CH3O2 = 46
  integer(ik4) , parameter :: ind_ISOPRD = 47
  integer(ik4) , parameter :: ind_ANO2 = 48
  integer(ik4) , parameter :: ind_ROOH = 49
  integer(ik4) , parameter :: ind_RO2 = 50
  integer(ik4) , parameter :: ind_ONIT = 51
  integer(ik4) , parameter :: ind_HO2 = 52
  integer(ik4) , parameter :: ind_O3 = 53
  integer(ik4) , parameter :: ind_OH = 54
  integer(ik4) , parameter :: ind_NO = 55
  integer(ik4) , parameter :: ind_NO2 = 56
  integer(ik4) , parameter :: ind_NO3 = 57
  integer(ik4) , parameter :: ind_C2O3 = 58

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
