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
!      write(stdout,*) '         _________________________________ '
!      write(stdout,*) '        (                                 ('     
!      write(stdout,*) '        )  RRRR  RRRRR  RRR   RRR  R   R  )'
!      write(stdout,*) '        (  R   R R     R   R R   R RR RR  ('
!      write(stdout,*) '        )  R   R R     R     R     R R R  )'
!      write(stdout,*) '        (  RRRR  RRRR  R RRR R     R R R  ('
!      write(stdout,*) '        )  R R   R     R   R R     R   R  )'
!      write(stdout,*) '        (  R  R  R     R   R R   R R   R  ('
!      write(stdout,*) '        )  R   R RRRRR  RRR   RRR  R   R  )'
!      write(stdout,*) '        (                                 ('
!      write(stdout,*) '        )         YOU DID IT,             )'
!      write(stdout,*) '        (         CONGRATULATIONS!        ('
!      write(stdout,*) '        )                                 )'
!      write(stdout,*) '        (   OOO   OOO   OOO  OOOOO  OOO   ('
!      write(stdout,*) '        )  O   O O   O O   O   O   O   O  )'
!      write(stdout,*) '        (  O   O O   O O       O   O      ('
!      write(stdout,*) '        )  O   O OOOOO  OOO    O    OOO   )'
!      write(stdout,*) '        (  O   O O   O     O   O       O  ('
!      write(stdout,*) '        )  O   O O   O O   O   O   O   O  )'
!      write(stdout,*) '        (   OOO  O   O  OOO  OOOOO  OOO   ('
!      write(stdout,*) '        )_________________________________)'
      write(stdout,*) ''
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
    write(stdout,*) '     [ Quentin Desmet,                 ]'
    write(stdout,*) '     [ quentin.desmet@legos.obs-mip.fr ]'
    write(stdout,*) '       Last update of the interface:'
    write(stdout,*) '       2022-02-10'
    write(stdout,*) ''
  end subroutine oasisxregcm_contacts

end module mod_oasis_signature
!
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
