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
module mod_mklaisai
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mklaisai

  character(len=16) , parameter :: timedim = 'time'
  character(len=16) , parameter :: varname1 = 'MONTHLY_LAI'
  character(len=16) , parameter :: varname2 = 'MONTHLY_SAI'
  character(len=32) , parameter :: varname3 = 'MONTHLY_HEIGHT_BOT'
  character(len=32) , parameter :: varname4 = 'MONTHLY_HEIGHT_TOP'

  contains

  subroutine mklaisai(laisaifile, &
                  monthly_lai,monthly_sai,monthly_top,monthly_bot)
    implicit none
    character(len=*) , intent(in) :: laisaifile
    real(rkx) , dimension(:,:,:,:) , intent(out) :: monthly_sai , monthly_lai
    real(rkx) , dimension(:,:,:,:) , intent(out) :: monthly_top , monthly_bot
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//laisaifile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,i_band)
    call gfread(gfile,varname1,monthly_lai,h_missing_value)
    call gfread(gfile,varname2,monthly_sai,h_missing_value)
    call gfread(gfile,varname3,monthly_top,h_missing_value)
    call gfread(gfile,varname4,monthly_bot,h_missing_value)
    call gfclose(gfile)
  end subroutine mklaisai

end module mod_mklaisai
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
