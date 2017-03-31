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
module mod_mkpopd
#ifdef CN
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkpopd , mkpopd_init , mkpopd_close

  character(len=16) , parameter :: varname = 'PDENS'

  type(globalfile) :: gfile

  contains

  subroutine mkpopd_init(popdfile)
    implicit none
    character(len=*) , intent(in) :: popdfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//popdfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,i_band)
  end subroutine mkpopd_init

  subroutine mkpopd(popd,it)
    implicit none
    real(rkx) , dimension(:,:) , intent(out) :: popd
    integer(ik4) , intent(in) :: it
    call gfread(gfile,varname,popd,it,0.1_rkx)
  end subroutine mkpopd

  subroutine mkpopd_close
    implicit none
    call gfclose(gfile)
  end subroutine mkpopd_close
#endif

end module mod_mkpopd
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
