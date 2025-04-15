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

#ifdef OASIS

module mod_oasis_signature

  use mod_mppparam
  use mod_stdio
  use mod_dynparam

  implicit none

  private

  public :: oasisxregcm_header , oasisxregcm_endscreen

  contains

  subroutine oasisxregcm_header
    implicit none
    !--------------------------------------------------------------------------
    if ( myid == italk ) then
      write(stdout,*) ''
      write(stdout,*) ' The code has been compiled using -DOASIS flag'
      write(stdout,*) '  (--enable-oasis with the program configure),'
      write(stdout,*) '  which implies the use of OASIS coupler.'
      write(stdout,*) '  Documentation should be available in RegCM/Doc/.'
      write(stdout,*) ''
      call oasisxregcm_contacts
      write(stdout,*) ''
    end if
  end subroutine oasisxregcm_header

  subroutine oasisxregcm_endscreen
    implicit none
    !--------------------------------------------------------------------------
    if ( myid == italk ) then
      write(stdout,*) '         _________________________________ '
      write(stdout,*) '        (                                 ('
      write(stdout,*) '        )  RRRR  RRRRR  RRR   RRR  R   R  )'
      write(stdout,*) '        (  R   R R     R   R R   R RR RR  ('
      write(stdout,*) '        )  R   R R     R     R     R R R  )'
      write(stdout,*) '        (  RRRR  RRRR  R RRR R     R R R  ('
      write(stdout,*) '        )  R R   R     R   R R     R   R  )'
      write(stdout,*) '        (  R  R  R     R   R R   R R   R  ('
      write(stdout,*) '        )  R   R RRRRR  RRR   RRR  R   R  )'
      write(stdout,*) '        (                                 ('
      write(stdout,*) '        )         YOU DID IT,             )'
      write(stdout,*) '        (         CONGRATULATIONS!        ('
      write(stdout,*) '        )                                 )'
      write(stdout,*) '        (   OOO   OOO   OOO  OOOOO  OOO   ('
      write(stdout,*) '        )  O   O O   O O   O   O   O   O  )'
      write(stdout,*) '        (  O   O O   O O       O   O      ('
      write(stdout,*) '        )  O   O OOOOO  OOO    O    OOO   )'
      write(stdout,*) '        (  O   O O   O     O   O       O  ('
      write(stdout,*) '        )  O   O O   O O   O   O   O   O  )'
      write(stdout,*) '        (   OOO  O   O  OOO  OOOOO  OOO   ('
      write(stdout,*) '        )_________________________________)'
      write(stdout,*) ''
    end if
  end subroutine oasisxregcm_endscreen

  subroutine oasisxregcm_contacts
    implicit none
    !--------------------------------------------------------------------------
    write(stdout,*) ''
    write(stdout,*) '       For any question, bug report, etc.'
    write(stdout,*) '       related to the use of OASIS in RegCM,'
    write(stdout,*) '       please contact:'
    write(stdout,*) '     [ Quentin Desmet,              ]'
    write(stdout,*) '     [ quentin.desmet@univ-tlse3.fr ]'
    write(stdout,*) '       Last update of the interface:'
    write(stdout,*) '       2023-09-28'
    write(stdout,*) ''
  end subroutine oasisxregcm_contacts

end module mod_oasis_signature
!
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
