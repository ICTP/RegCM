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

module mod_ecearth_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date

  private

  public :: echvars , echcmorvars
  public :: find_ecearth_sst , find_ecearth_topo
  public :: find_ecearth_dim , find_ecearth_file

  integer(ik4) , parameter :: nvars = 6
  character(len=3) , target , dimension(nvars) :: echvars = &
            ['t  ' , 'z  ' , 'q  ' , 'u  ' , 'v  ', 'XXX']
  character(len=3) , target , dimension(nvars) :: echcmorvars = &
            ['ta ' , 'XXX' , 'hus' , 'ua ' , 'va ', 'aps']

  contains

  subroutine find_ecearth_sst(fname,idate,cmor)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    logical , intent(in) :: cmor
    character(len=10) :: d1 , d2
    integer(ik4) :: y , m , d , h
    integer(ik4) :: y1 , y2 , m1 , m2
    if ( cmor ) then
      call split_idate(idate,y,m,d,h)
      y1 = y
      m1 = m
      if ( d == 1 .and. h == 0 ) then
        m1 = m - 1
        if ( m1 == 0 ) then
          m1 = 12
          y1 = y1 - 1
        end if
      end if
      y2 = y1
      m2 = m1 + 1
      if ( m2 > 12 ) then
        m2 = 1
        y2 = y2 + 1
      end if
      write(d1,'(i0.4,i0.2,i0.2,i0.2)') y1, m1, 1, 6
      write(d2,'(i0.4,i0.2,i0.2,i0.2)') y2, m2, 1, 0
      if ( idate < 2005120100 ) then
        fname = trim(inpglob)//pthsep//'EC-EARTH'//pthsep//'SST'// &
                pthsep//'sst_6hrLev_ICHEC-EC-EARTH_historical'// &
                '_r12i1p1_'//d1//'00-'//d2//'00.nc'
      else
        fname = trim(inpglob)//pthsep//'EC-EARTH'//pthsep//'SST'// &
                pthsep//'sst_6hrLev_ICHEC-EC-EARTH_rcp'//ssttyp(4:5)//  &
                '_r12i1p1_'//d1//'00-'//d2//'00.nc'
      end if
    else
      if ( .not. date_in_scenario(idate,5,.true.) ) then
        fname = trim(inpglob)//'/EC-EARTH/SST/RF/ich1_sst_1950-2009.nc'
      else
        fname = trim(inpglob)//'/EC-EARTH/SST/RCP'//ssttyp(4:5)//&
                  '/ic'//ssttyp(4:4)//'1_sst_2006-2100.nc'
      end if
    end if
  end subroutine find_ecearth_sst

  subroutine find_ecearth_topo(fname,cmor)
    implicit none
    character(len=256) , intent(out) :: fname
    logical , intent(in) :: cmor
    if ( cmor ) then
      fname = trim(inpglob)// &
                '/EC-EARTH/fixed/orog_fx_EC-EARTH_historical_r0i0p0.nc'
    else
      fname = trim(inpglob)//'/EC-EARTH/fixed/ecearth.nc'
    end if
  end subroutine find_ecearth_topo

  subroutine find_ecearth_dim(dim_filename,cmor)
    implicit none
    character(len=256) , intent(out) :: dim_filename
    logical , intent(in) :: cmor
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    if ( cmor ) then
      dim_filename = trim(inpglob)// &
                '/EC-EARTH/RF/ta/ta_6hrLev_ICHEC-EC-EARTH_'//&
                'historical_r12i1p1_196901010600-196902010000.nc'
    else
      dim_filename = trim(inpglob)//'/EC-EARTH/fixed/ecearth.nc'
    end if
  end subroutine find_ecearth_dim

  subroutine find_ecearth_file(ecearth_filename,var,idate,cmor)
    implicit none
    logical , intent(in) :: cmor
    character(len=256) , intent(out) :: ecearth_filename
    character(len=*) , intent(in) :: var
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=256) :: inname
    character(len=10) :: d1 , d2
    integer(ik4) :: y , m , d , h
    integer(ik4) :: y1 , y2 , m1 , m2
    call split_idate(idate,y,m,d,h)
    if ( cmor ) then
      y1 = y
      m1 = m
      if ( d == 1 .and. h == 0 ) then
        m1 = m - 1
        if ( m1 == 0 ) then
          m1 = 12
          y1 = y1 - 1
        end if
      end if
      y2 = y1
      m2 = m1 + 1
      if ( m2 > 12 ) then
        m2 = 1
        y2 = y2 + 1
      end if
      write(d1,'(i0.4,i0.2,i0.2,i0.2)') y1, m1, 1, 6
      write(d2,'(i0.4,i0.2,i0.2,i0.2)') y2, m2, 1, 0
      if ( idate < 2005120100 ) then
        ecearth_filename = trim(inpglob)//pthsep//'EC-EARTH'//pthsep//'RF'// &
                pthsep//trim(var)//pthsep//trim(var)// &
                '_6hrLev_ICHEC-EC-EARTH_historical_r12i1p1_'// &
                d1//'00-'//d2//'00.nc'
      else
        ecearth_filename = trim(inpglob)//pthsep//'EC-EARTH'//pthsep// &
                'RCP'//dattyp(4:5)//pthsep//trim(var)//pthsep//trim(var)// &
                '_6hrLev_ICHEC-EC-EARTH_historical_r12i1p1_'// &
                d1//'00-'//d2//'00.nc'
      end if
    else
      if ( m == 1 .and. d == 1 .and. h == 0 ) y = y - 1
      if ( .not. date_in_scenario(idate,5,.true.) ) then
        write (inname,'(a,a,i0.4,a,a,a,i0.4,a)') 'RF', pthsep, y, pthsep, &
                   'ich1_', trim(var)//'_', y, '.nc'
      else
        write (inname,'(a,a,i0.4,a,a,a,i0.4,a)') 'RCP'//dattyp(4:5),  &
                   pthsep, y, pthsep, 'ich1_', trim(var)//'_', y, '.nc'
      end if
      ecearth_filename = trim(inpglob)//'/EC-EARTH/'//inname
    end if
  end subroutine find_ecearth_file

end module mod_ecearth_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
