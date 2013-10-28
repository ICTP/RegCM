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
!
!
program lakeprog
  use mod_lake
  implicit none
!
  real (8) , parameter :: dpth = 50.0D0 ! Lake depth in m
  real (8) , parameter :: dt = 300.0D0   ! Lake internal timestep in sec
  real (8) , parameter :: xl = 45.0D0    ! Latitude of point in degrees
  real (8) , parameter :: zl = 10.0D0    ! Elevation at which vl is computed
!
! First time into, try to reach equilibrium calling lake this number of times
!
  integer , parameter :: nconv = 200
!
  integer , parameter :: inunit = 111
  integer , parameter :: iounit = 222
  integer :: iloop , idp , iconv
  character(256) :: inpfile , outfile , prgname

  real (8) :: tl , vl , ql , fsw , flw , hsen , prec , evl

  call get_command_argument(0,value=prgname)
  call get_command_argument(1,value=inpfile)
  call get_command_argument(2,value=outfile)

  open (inunit,file=inpfile,status='old',action='read', &
        form='formatted',err=100)
  open (iounit,file=outfile,status='replace',action='write', &
        form='formatted',err=200)

  call initlake(dpth)

  ! Skip header line in input file
  read(inunit, *)
  iloop = 0
  read(inunit, '(8f9.3)', end=999) tl, vl, ql, fsw, flw, hsen, prec, evl
  evl = evl*1.0D-3

  ! Try to find an equilibrium
  do iconv = 1 , nconv
    call lake(dt,tl,vl,zl,ql,fsw,flw,hsen,xl,prec,evl)
  end do
  call writeout

  ! Simulate now
  iloop = 1
  do
    read(inunit, '(8f9.3)', end=999) tl, vl, ql, fsw, flw, hsen, prec, evl
    evl = evl*1.0D-3
    call lake(dt,tl,vl,zl,ql,fsw,flw,hsen,xl,prec,evl)
    call writeout
  end do

999 continue
  stop 'Successful end'

100 continue
  write(6, *) 'Cannot open input file ', trim(inpfile)
  call usage
200 continue
  write(6, *) 'Cannot open output file ', trim(outfile)
  call usage

  contains

  subroutine usage
    implicit none
    write (6, *) trim(prgname), ' : 1D lake model'
    write (6, *) 'Usage : '
    write (6, *) '      ',trim(prgname), ' inputfile outputfile'
    write (6, *) ' '
    write (6, *) 'Input format is 8f9.3'
    write (6, *) 'Input is: tl, vl, ql, fsw, flw, hsen, prec, evl'
    write (6, *) 'tl   = air temperature in K'
    write (6, *) 'vl   = wind speed in m/s'
    write (6, *) 'ql   = humidity mixing ratio in kg/kg'
    write (6, *) 'fsw  = short wave soil flux W/m^2'
    write (6, *) 'flw  = long wave soil flux W/m^2'
    write (6, *) 'hsen = sensible heat flux W/m^2'
    write (6, *) 'prec = precipitation in mm/sec'
    write (6, *) 'evl  = evaporation in mm/sec'

    stop 'Invalid input'
  end subroutine usage

  subroutine writeout
    implicit none
    write(iounit, '(a)') '###########################################'
    write(iounit, '(a,i10)')     'Loop index  : ', iloop
    write(iounit, '(a,f12.4,a)') 'Ground Temp : ', tgl, 'K'
    write(iounit, '(a,f12.4,a)') 'Ice depth   : ', aveice, 'mm'
    write(iounit, '(a,f12.4,a)') 'Snow depth  : ', hsnow, 'mm'
    write(iounit, '(a,f12.4,a)') 'Evaporation : ', evl, 'mm/sec'
    write(iounit, '(a,i3,a)')    'Temperature profile from 1 to ',idep,'m'
    do idp = 1 , idep
      write(iounit, '(i3,f12.4)') idp, tprof(idp)
    end do
    write(iounit, '(a)') '###########################################'
  end subroutine writeout

end program lakeprog
