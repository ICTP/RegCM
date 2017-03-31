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
module mod_mkvocef
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkvocef

  integer , parameter :: nvocs = 6

  character(len=16) , parameter , dimension(nvocs):: varname = &
          (/ 'ef_btr' , 'ef_crp' , 'ef_fdt' , &
             'ef_fet' , 'ef_grs' , 'ef_shr'/)

  contains

  subroutine mkvocef(vocfile,vocef)
    implicit none
    character(len=*) , intent(in) :: vocfile
    real(rkx) , dimension(:,:,:) , intent(out) :: vocef
    integer(ik4) :: n
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//vocfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,i_band)
    do n = 1 , nvocs
      call gfread(gfile,varname(n),vocef(:,:,n),h_missing_value)
    end do
    call gfclose(gfile)
  end subroutine mkvocef

end module mod_mkvocef
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
