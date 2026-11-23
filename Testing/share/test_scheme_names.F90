!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Unit test of Share/mod_scheme_names.F90: the codes of the scheme names
! of every kind (the integer codes of the model), case and separator
! insensitivity, aliases, the removed schemes, unknown names, and the
! name of each code back.

program test_scheme_names
  use mod_intkinds, only : ik4
  use mod_scheme_names, only : scheme_code, scheme_name, scheme_list, &
    unknown_scheme
  use share_test_utils, only : check, test_end
  implicit none
  character(len=256) :: msg

  call check(scheme_code('pbl','holtslag',msg) == 1,'pbl holtslag')
  call check(scheme_code('pbl','uw',msg) == 2,'pbl uw')
  call check(scheme_code('pbl','gfs',msg) == 3,'pbl gfs')
  call check(scheme_code('pbl','myj',msg) == 4,'pbl myj')
  call check(scheme_code('pbl','shin-hong',msg) == 5,'pbl shin-hong')
  call check(scheme_code('pbl','none',msg) == 0,'pbl none')
  call check(scheme_code('cumulus','grell',msg) == 2,'cumulus grell')
  call check(scheme_code('cumulus','emanuel',msg) == 4,'cumulus emanuel')
  call check(scheme_code('cumulus','tiedtke',msg) == 5,'cumulus tiedtke')
  call check(scheme_code('cumulus','kain-fritsch',msg) == 6, &
    'cumulus kain-fritsch')
  call check(scheme_code('cumulus','shallow',msg) == -1,'cumulus shallow')
  call check(scheme_code('microphysics','noto',msg) == 2,'micro noto')
  call check(scheme_code('microphysics','wsm5',msg) == 3,'micro wsm5')
  call check(scheme_code('microphysics','wsm7',msg) == 4,'micro wsm7')
  call check(scheme_code('microphysics','wdm7',msg) == 5,'micro wdm7')
  call check(scheme_code('cloud_fraction','sundqvist',msg) == 0, &
    'cloud sundqvist')
  call check(scheme_code('cloud_fraction','xu-randall',msg) == 1, &
    'cloud xu-randall')
  call check(scheme_code('cloud_fraction','digital',msg) == 7, &
    'cloud digital')
  call check(scheme_code('radiation','ccm',msg) == 0,'radiation ccm')
  call check(scheme_code('radiation','rrtmg',msg) == 1,'radiation rrtmg')
  ! case, separators, blanks and aliases
  call check(scheme_code('cumulus','  Kain_Fritsch',msg) == 6, &
    'case and separator insensitive')
  call check(scheme_code('cumulus','KF',msg) == 6,'alias kf')
  call check(scheme_code('pbl','ShinHong',msg) == 5,'alias shinhong')
  call check(scheme_code('cloud_fraction','xuran',msg) == 1,'alias xuran')
  ! removed and unknown
  call check(scheme_code('cumulus','kuo',msg) == unknown_scheme .and. &
    index(msg,'historical') > 0,'kuo removed')
  call check(scheme_code('cumulus','betts-miller',msg) == unknown_scheme &
    .and. index(msg,'historical') > 0,'betts-miller removed')
  call check(scheme_code('microphysics','subex',msg) == unknown_scheme &
    .and. index(msg,'historical') > 0,'subex microphysics removed')
  call check(scheme_code('cloud_fraction','subex',msg) == unknown_scheme &
    .and. index(msg,'sundqvist') > 0,'subex cloud fraction renamed')
  call check(scheme_code('pbl','tiedtke',msg) == unknown_scheme .and. &
    index(msg,'holtslag') > 0,'unknown name lists the schemes')
  call check(scheme_code('turbulence','uw',msg) == unknown_scheme, &
    'unknown kind')
  ! back from the codes
  call check(scheme_name('cumulus',6) == 'kain-fritsch','name of cumulus 6')
  call check(scheme_name('microphysics',2) == 'noto','name of micro 2')
  call check(scheme_name('cloud_fraction',0) == 'sundqvist', &
    'name of cloud fraction 0')
  call check(scheme_name('cumulus',1) == 'unknown','no name for code 1')
  call check(index(scheme_list('radiation'),'ccm, rrtmg') == 1, &
    'list of radiation schemes')
  call test_end('test_scheme_names')
end program test_scheme_names
