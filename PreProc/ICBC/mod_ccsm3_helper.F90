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

module mod_ccsm3_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_message

  private

  public :: find_ccsm3_topo , find_ccsm3_file

  character(len=32) , parameter :: ccsm3_experiment = 'OCN.199904.E30.cam2.h1.'
  character(len=64) , parameter :: ccsm3_init = &
          'OCN.199904.E30.cam2.i.1999-05-01-00000.nc'
  character(len=16) :: ccsm3_date

  contains

  subroutine find_ccsm3_topo(topo_filename)
    implicit none
    character(len=256) , intent(out) :: topo_filename
    topo_filename = trim(inpglob)//'/CCSM3/'//trim(ccsm3_init)
  end subroutine find_ccsm3_topo

  subroutine find_ccsm3_file(ccsm3_filename,y,m,d,h)
    implicit none
    character(len=256) , intent(out) :: ccsm3_filename
    integer(ik4) , intent(in) :: y,m,d,h
    integer(ik4)  :: yy , mm , dd , hh , icheck , inow , ii
    integer(ik4)  , parameter :: oneyear = 10000000
    logical :: iexist
    yy = y
    mm = m
    dd = d
    hh = h
    ! Lucky search
    write(ccsm3_date,'(i0.4,a1,i0.2,a1,i0.2,a1,i0.5)') &
           yy , '-', mm, '-', dd, '-', hh*3600
    ccsm3_filename = trim(inpglob)//'/CCSM3/'//trim(ccsm3_experiment)// &
               ccsm3_date//'.nc'
    inquire(file=ccsm3_filename,exist=iexist)
    if ( iexist ) return
    ! Search if file for today is here
    do ii = hh , 0 , -1
      write(ccsm3_date,'(i0.4,a1,i0.2,a1,i0.2,a1,i0.5)') &
             yy , '-', mm, '-', dd, '-', ii*3600
      ccsm3_filename = trim(inpglob)//'/CCSM3/'//trim(ccsm3_experiment)// &
                 ccsm3_date//'.nc'
      inquire(file=ccsm3_filename,exist=iexist)
      if ( iexist ) return
    end do
    call getbackoneday(yy,mm,dd)
    ! We will search data file from the requested date back for 2 years.
    icheck = y*1000000+m*10000+d*100+h
    do
      inow = yy*1000000+mm*10000+dd*100+hh
      do ii = 23 , 0 , -1
        write(ccsm3_date,'(i0.4,a1,i0.2,a1,i0.2,a1,i0.5)') &
           yy , '-', mm, '-', dd, '-', ii*3600
        ccsm3_filename = trim(inpglob)//'/CCSM3/'//trim(ccsm3_experiment)// &
               ccsm3_date//'.nc'
        inquire(file=ccsm3_filename,exist=iexist)
        if ( iexist ) return
      end do
      if (icheck-oneyear > inow ) exit
      call getbackoneday(yy,mm,dd)
    end do
    call die('gnhnc_sst','ccsm3 file with requested data cannot be found',1)
  end subroutine find_ccsm3_file
!
  subroutine getbackoneday(y,m,d)
    implicit none
    integer(ik4) , intent(inout) :: y , m , d
    integer(ik4) , dimension(12) , parameter :: dpm = &
            (/31,28,31,30,31,30,31,31,30,31,30,31/)
    d = d - 1
    if ( d == 0 ) then
      m = m - 1
      d = dpm(m)
      if ( m == 0 ) then
        m = 12
        d = 31
        y = y -1
      end if
    end if
  end subroutine getbackoneday

end module mod_ccsm3_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
