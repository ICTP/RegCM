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

module mod_hadgem_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date

  private

  public :: havars , habase1 , habase2 , habase3
  public :: find_hadgem_sst
  public :: find_hadgem_dim , find_hadgem_file
  public :: find_hadgem_ufile , find_hadgem_vfile

  character(len=16) :: hadgem_date

  integer(ik4) , parameter :: nvars = 6
  character(len=3) , target , dimension(nvars) :: havars = &
            (/'ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps '/)

  character(len=64) :: habase1  = '_6hrLev_HadGEM2-ES_historical'
  character(len=64) :: habase2  = '_6hrLev_HadGEM2-ES_rcp'
  character(len=64) :: habase3 = '_r1i1p1_'

  contains

  subroutine find_hadgem_sst(fname,idate)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    if ( .not. date_in_scenario(idate,5) ) then
      if ( idate > 1959110100 ) then
        fname = trim(inpglob)//pthsep//'HadGEM2'//pthsep//'SST'// &
                pthsep//'tos_Omon_HadGEM2-ES_historical'// &
                '_r1i1p1_195912-200512.nc'
      else
        fname = trim(inpglob)//pthsep//'HadGEM2'//pthsep//'SST'// &
                pthsep//'tos_Omon_HadGEM2-ES_historical'// &
                '_r1i1p1_185912-195911.nc'
      end if
    else
      if ( idate < 2099110100 ) then
        fname = trim(inpglob)//pthsep//'HadGEM2'//pthsep//'SST'// &
                pthsep//'tos_Omon_HadGEM2-ES_rcp'//ssttyp(4:5)//  &
                '_r1i1p1_200512-209911.nc'
      else
        fname = trim(inpglob)//pthsep//'HadGEM2'//pthsep//'SST'// &
                pthsep//'tos_Omon_HadGEM2-ES_rcp'//ssttyp(4:5)//  &
                '_r1i1p1_209912-219911.nc'
      end if
    end if
  end subroutine find_hadgem_sst

  subroutine assemble_path(fname,scen,var,d1,d2)
    implicit none
    character(len=256) , intent(out) :: fname
    character(len=*) , intent(in) :: scen
    character(len=*) , intent(in) :: var
    character(len=*) , intent(in) :: d1
    character(len=*) , intent(in) :: d2
    if ( scen == 'RF' ) then
      fname = trim(inpglob)//pthsep//'HadGEM2'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(habase1)//  &
              trim(habase3)//trim(d1)//'-'//trim(d2)//'.nc'
    else
      fname = trim(inpglob)//pthsep//'HadGEM2'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(habase2)//  &
              scen(4:5)//trim(habase3)//trim(d1)//'-'//trim(d2)//'.nc'
    end if
  end subroutine assemble_path

  subroutine find_hadgem_dim(dim_filename)
    implicit none
    character(len=256) , intent(out) :: dim_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(dim_filename,'RF','ta','1990120106','1991030100')
  end subroutine find_hadgem_dim

  subroutine find_hadgem_ufile(ufile_filename)
    implicit none
    character(len=256) , intent(out) :: ufile_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(ufile_filename,'RF','ua','1990120106','1991030100')
  end subroutine find_hadgem_ufile

  subroutine find_hadgem_vfile(vfile_filename)
    implicit none
    character(len=256) , intent(out) :: vfile_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(vfile_filename,'RF','va','1990120106','1991030100')
  end subroutine find_hadgem_vfile

  subroutine find_hadgem_file(hadgem_filename,var,idate)
    implicit none
    character(len=256) , intent(out) :: hadgem_filename
    character(len=*) , intent(in) :: var
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=10) :: d1 , d2
    integer(ik4) :: y , m , d , h
    integer(ik4) :: yy , mm
    integer(ik4) :: icheck , inow
    call split_idate(idate,y,m,d,h)
    select case (var)
      case ('ps')
        inow = y*1000000+m*10000+d*100+h
        if ( .not. date_in_scenario(idate,5,.false.) .or. &
             inow == 2005120100 ) then
          icheck = y*1000000+120200
          if ( inow > icheck ) y = y + 1
          write(d1,'(i0.4i0.2i0.2i0.2)') y-1, 12, 2, 6
          write(d2,'(i0.4i0.2i0.2i0.2)') y, 12, 2, 0
          call assemble_path(hadgem_filename,'RF',havars(6),d1,d2)
        else
          icheck = y*1000000+120100
          if ( inow > icheck ) y = y + 1
          write(d1,'(i0.4i0.2i0.2i0.2)') y-1, 12, 1, 6
          write(d2,'(i0.4i0.2i0.2i0.2)') y, 12, 1, 0
          call assemble_path(hadgem_filename,'RCP'//dattyp(4:5), &
                             havars(6),d1,d2)
        end if
      case default
        if ( dattyp(4:5) == '26' ) then
          inow = y*1000000+m*10000+d*100+h
          if ( .not. date_in_scenario(idate,5,.false.) .or. &
               inow == 2005120100 ) then
            yy = y
            mm = (m/3)*3
            if ( mm == 0 ) then
              yy = y - 1
              mm = 12
            end if
            icheck = yy*1000000+mm*10000+106
            if ( inow < icheck ) then
              mm = mm - 3
              if ( mm == 0 ) then
                yy = yy - 1
                mm = 12
              end if
            end if
            write(d1,'(i0.4i0.2i0.2i0.2)') yy, mm, 1, 6
            mm = mm + 3
            if ( mm > 12 ) then
              mm = 3
              yy = yy + 1
            end if
            write(d2,'(i0.4i0.2i0.2i0.2)') yy, mm, 1, 0
            call assemble_path(hadgem_filename,'RF',var,d1,d2)
          else
            icheck = y*1000000+120100
            if ( inow > icheck ) y = y + 1
            write(d1,'(i0.4i0.2i0.2i0.2)') y-1, 12, 1, 6
            write(d2,'(i0.4i0.2i0.2i0.2)') y, 12, 1, 0
            call assemble_path(hadgem_filename,'RCP'//dattyp(4:5),var,d1,d2)
          end if
        else
          inow = y*1000000+m*10000+d*100+h
          yy = y
          mm = (m/3)*3
          if ( mm == 0 ) then
            yy = y - 1
            mm = 12
          end if
          if ( .not. date_in_scenario(idate,5,.false.) .or. &
               inow == 2005120100 ) then
            icheck = yy*1000000+mm*10000+106
            if ( inow < icheck ) then
              mm = mm - 3
              if ( mm == 0 ) then
                yy = yy - 1
                mm = 12
              end if
            end if
            write(d1,'(i0.4i0.2i0.2i0.2)') yy, mm, 1, 6
            mm = mm + 3
            if ( mm > 12 ) then
              mm = 3
              yy = yy + 1
            end if
            write(d2,'(i0.4i0.2i0.2i0.2)') yy, mm, 1, 0
            call assemble_path(hadgem_filename,'RF',var,d1,d2)
          else
            icheck = yy*1000000+mm*10000+106
            if ( inow < icheck ) then
              mm = mm - 3
              if ( mm == 0 ) then
                yy = yy - 1
                mm = 12
              end if
            end if
            write(d1,'(i0.4i0.2i0.2i0.2)') yy, mm, 1, 6
            mm = mm + 3
            if ( mm > 12 ) then
              mm = 3
              yy = yy + 1
            end if
            write(d2,'(i0.4i0.2i0.2i0.2)') yy, mm, 1, 0
            call assemble_path(hadgem_filename,'RCP'//dattyp(4:5),var,d1,d2)
          end if
        end if
    end select
  end subroutine find_hadgem_file

end module mod_hadgem_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
