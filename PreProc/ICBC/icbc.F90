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

#ifdef PNETCDF
subroutine myabort
  use mod_stdio
  use mpi
  implicit none
  integer :: ierr
  write(stderr,*) ' Execution terminated because of runtime error'
  call mpi_abort(mpi_comm_self,1,ierr)
end subroutine myabort
#else
subroutine myabort
  implicit none
  stop ' Execution terminated because of runtime error'
end subroutine myabort
#endif

program icbc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                      !
!   ICBC is the third component of the REGional Climate Modeling       !
!   (RegCM) system version 3.0 and used to access archived global      !
!   analysed datasets at regular latitude-longititude (NNRP1, NNRP2,   !
!   ERA40, ERAIN,EIN75,EIN15, EIN25, GFS11)                            !
!   or original T159 (N80) datasets(ERAHI),                            !
!   or T42 datasets at Gaussian grids (ECWCRP, simply ECMWF), as well  !
!   as NEST run from previous FVGCM run (FVGCM), ECHAM5 run  (EH5OM)   !
!   and RegCM run (FNEST).                                             !
!                                                                      !
!   The present ICBC code could treat NNRP1, NNRP2, ECWCRP, ERA40,     !
!   ERAIN, EIN75, EIN15, EIN25, GFS11, ERAHI, FVGCM, EH5OM, ECEXY,     !
!   and RegCM datasets,  4 times daily.                                !
!                                                                      !
!                        Xunqiang Bi, ESP group, Abdus Salam ICTP      !
!                                                October 07, 2009      !
!                                                                      !
!   CCSMN: unpacked CCSM NETCDF (six hourly) data                      !
!          Get from EarthSystemGrid                                    !
!   CAM4N: CAM2/4 NETCDF (six hourly) data                             !
!   CA_XX: Canadian model CMIP5 netCDF (6 hourly) data                 !
!          The XX stands for one in RF, 26, 45, 60, 85                 !
!          which reference the RCP scenarios                           !
!   HA_XX: The Hadgem model CMIP5 netCDF (6 hourly) data               !
!          The XX stands for one in RF, 26, 45, 60, 85                 !
!          which reference the RCP scenarios                           !
!   E_ICH: The EC-Earth dataset (6 hourly) on pressure levels          !
!          The original grib on sperical coordinates was processed     !
!          with the script to be found in the Tools/Script directory   !
!   IP_XX: The IPSL CMIP5 dataset (6 hourly)                           !
!          The XX stands for one in RF, 26, 45, 60, 85                 !
!          which reference the RCP scenarios                           !
!   NNRP1: NCEP/NCAR Reanalysis datasets are available at:             !
!          ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis/            !
!          Current holdings: 1948 - present, 2.5x2.5L13, netCDF.       !
!   NNRP2: NCEP/DOE AMIP-II Reanalysis (Reanalysis-2) are at:          !
!          ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis2/           !
!          Current holdings: 1979 - 2009, 2.5x2.5L13, netCDF.          !
!   CFSXX: CFS seasonal forecast on pressure level.                    !
!          XX stands for ensemble member (in 01 02 03 04)              !
!   ECMWF: ECMWF TOGA/WCRP Uninitialized Data - (ECWCRP)               !
!          NCAR MSS:/TRENBERT/CTEC/ , ET42yymmdd, where yy = year,     !
!          mm = month, dd = day = 01,04,07,10,13,16,19,22,25,28, or 31 !
!          Current holdings: January, 1993 - December, 1997            !
!          Reformatted by PWC/ICTP to direct-access binary,            !
!          T42L15, Gaussian Grid.                                      !
!   EH5OM: EH5OM run by the MPI at Hamburg, T63, Gaussian grid.        !
!          For present day  run: 1941 - 2000;                          !
!          For A1B scenario run: 2001 - 2100.                          !
!          17 pressure levels, 4 times daily, direct-access binary.    !
!   ERA40: ECMWF 40 year reanalysis datasets are available at:         !
!          http://data.ecmwf.int/data/d/era40_daily/                   !
!          Current holdings: 01/09/1957 - 31/08/2002,                  !
!          Pressure levels, 2.5x2.5L23, 4 times daily.                 !
!   ERAIN/EIN15: ECMWF INTERIM 10 year reanalysis datasets             !
!          Current holdings: 01/01/1989 - 31/05/2009,                  !
!          Pressure levels, 1.5x1.5L37, 4 times daily.                 !
!   EIN75: ECMWF INTERIM 10 year reanalysis datasets                   !
!          Current holdings: 01/01/1989 - 31/12/2007,                  !
!          Pressure levels, 0.75x0.75L37, 4 times daily.               !
!   EIN25: ECMWF INTERIM 10 year reanalysis datasets                   !
!          Current holdings: 01/01/1989 - 31/12/1998,                  !
!          Pressure levels, 2.5x2.5L37, 4 times daily.                 !
!   GFS11: NCEP Global Forecast System (GFS) product FNL are           !
!                                                available at:         !
!          http://dss.ucar.edu/datasets/ds083.2/data/fnl-yyyymm/       !
!          Current holdings: 01/01/2000 - present,                     !
!          Pressure levels, 1.0x1.0L27, 4 times daily.                 !
!   ERAHI: ECMWF 40 year reanalysis datasets, origigal model level     !
!          fields: T, U, V and log(Ps) are in spectral coefficients    !
!          Oro and Q are at the reduced Gaussian grids.                !
!          T159L60 (N80L60), 01/09/1957 - 31/08/2002.                  !
!   ECEXY: ECMWF Ensemble forecast model. The X stands for model       !
!          version, the Y for the ensemble member number.              !
!          For example, ECE21 stands for ECMWF Ensemble model 2,       !
!          member number 1.                                            !
!   FVGCM: FVGCM run by the PWC group of Abdus Salam ICTP.             !
!          For present day run: 1960 - 1990;                           !
!          For A2          run: 2070 - 2100.                           !
!          1x1.25L18, 4 times daily, direct-access binary.             !
!   JRA55: The Japan Meteorological Agency (JMA) JRA-55, the second    !
!          Japanese global atmospheric reanalysis project.             !
!   FNEST: do Further oneway NESTing from previous RegCM run.          !
!                                                                      !
!   The code need NetCDF library to access ERA40, ERAIN (EIN75/15/25), !
!   NNRP1 and NNRP2 data.                                              !
!   And we have already provided the NetCDF libraries for many         !
!   platforms, if your platform is unfortunately out of our support,   !
!   you need install the netcdf library by yourself.                   !
!   The code need EMOSLIB library to access ERAHI (T159L60) data, we   !
!   have just provided EMOSLIB library for LINUX PGI5 and IBM AIX.     !
!                                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_message
  use mod_header
  use mod_stdio
  use mod_memutil
  use mod_mksst
  use mod_date
  use mod_grid
  use mod_date
  use mod_ecwcp
  use mod_eh5om
  use mod_ein
  use mod_era40
  use mod_era5
  use mod_erahi
  use mod_ecens
  use mod_fvgcm
  use mod_ncep
  use mod_nest
  use mod_gn6hnc
  use mod_write
  use mod_projections
#ifdef PNETCDF
  use mpi
#endif

  implicit none

  integer(ik4) :: nnn
  type(rcm_time_and_date) :: idate , iodate
  type(rcm_time_interval) :: tdiff , tbdy
  integer(ik4) :: nsteps
  integer(ik4) :: ierr
  character(len=256) :: namelistfile, prgname
  type(anyprojparams) :: pjpara

#ifdef PNETCDF
  call mpi_init(ierr)
#endif

  call header('icbc')
  !
  ! Read input global namelist
  !
  call get_command_argument(0,value=prgname)
  call get_command_argument(1,value=namelistfile)
  call initparam(namelistfile, ierr)
  if ( idynamic == 2 ) then
    write(stdout, *) 'Using non hydrostatic parameters'
    write(stdout, '(a,f10.2)') ' base_state_pressure    = ', base_state_pressure
    write(stdout, '(a,f10.2)') ' logp_lrate             = ', logp_lrate
  end if

  if ( idynamic == 3 ) then
    write(stdout, *) 'Using Moloch non-hydrostatic dynamical core'
  end if

  if ( dattyp == 'FVGCM' .or. dattyp == 'EH5RF' .or. &
       dattyp == 'EH5A2' .or. dattyp == 'EH5B1' .or. dattyp == 'EHA1B') then
    call init_globwindow(namelistfile,lat0,lon0,lat1,lon1)
  else if ( dattyp == 'FNEST' ) then
    call init_fnestparam(namelistfile,coarsedir,coarsedom)
  end if

  if ( ierr /= 0 ) then
    write ( stderr, * ) 'Parameter initialization not completed'
    write ( stderr, * ) 'Usage : '
    write ( stderr, * ) '          ', trim(prgname), ' regcm.in'
    write ( stderr, * ) ' '
    write ( stderr, * ) 'Check argument and namelist syntax'
    call die('icbc','Check argument and namelist syntax',1)
  end if
!
  call memory_init

  call init_grid(jx,iy,kz)
  call init_output

  pjpara%pcode = iproj
  pjpara%ds = ds*1000.0_rk8
  pjpara%clat = clat
  pjpara%clon = clon
  pjpara%plat = plat
  pjpara%plon = plon
  pjpara%trlat1 = truelatl
  pjpara%trlat2 = truelath
  pjpara%nlon = jx
  pjpara%nlat = iy
  pjpara%rotparam = .true.

  if ( idynamic == 3 ) then
    pjpara%staggerx = .true.
    pjpara%staggery = .false.
    call pju%initialize(pjpara)
    pjpara%staggerx = .false.
    pjpara%staggery = .true.
    call pjv%initialize(pjpara)
  else
    pjpara%staggerx = .true.
    pjpara%staggery = .true.
    call pjd%initialize(pjpara)
  end if

  if (dattyp == 'CCSMN' .or. dattyp == 'CAM4N' .or. dattyp == 'CCSM3' .or. &
      dattyp(1:3) == 'CA_' .or. dattyp(1:3) == 'GF_' .or. &
      dattyp(1:3) == 'NO_' .or. dattyp(1:3) == 'CC_' ) then
    if (ical /= noleap ) then
      write(stderr,*) 'Calendar should be set to noleap'
      call die('icbc','Calendar mismatch',1)
    end if
  end if
  if (dattyp(1:3) == 'HA_' ) then
    if ( ical /= y360 ) then
      write(stderr,*) 'Calendar should be set to 360_day'
      call die('icbc','Calendar mismatch',1)
    end if
  end if
  if (dattyp(1:2) == 'MP' .or. dattyp(1:2) == 'E5' .or. &
      dattyp(1:3) == 'LGM' ) then
    if ( ical /= gregorian ) then
      write(stderr,*) 'Calendar should be set to gregorian'
      call die('icbc','Calendar mismatch',1)
    end if
  end if

  tdiff = globidate2-globidate1
  tbdy = rcm_time_interval(ibdyfrq,uhrs)
  nsteps = nint(tohours(tdiff))/ibdyfrq + 1

  write (stdout,*) 'GLOBIDATE1 : ' , tochar(globidate1)
  write (stdout,*) 'GLOBIDATE2 : ' , tochar(globidate2)
  write (stdout,*) 'NSTEPS     : ' , nsteps

  idate = globidate1
  iodate = idate

  if ( dattyp(1:4) == 'NNRP' .or. dattyp(1:3) == 'CFS' ) then
    call init_ncep
  else if ( dattyp == 'ECMWF' ) then
    call init_ecwcp
  else if ( dattyp == 'ERA40' ) then
    call init_era40
  else if ( dattyp(1:4) == 'ERA5' ) then
    call init_era5
  else if ( dattyp == 'ERAIN' .or. dattyp(1:3) == 'EIN' .or. &
            dattyp == 'EIXXX' ) then
    call init_ein
  else if ( dattyp(1:3) == 'ECE' ) then
    call init_ecens
  else if ( dattyp == 'ERAHI' ) then
    call init_ehi
  else if ( dattyp(1:2) == 'EH' ) then
    call init_eh5om
  else if ( dattyp == 'FVGCM' ) then
    call init_fvgcm
  else if ( dattyp == 'FNEST' ) then
    call init_nest
  else
    if ( dattyp(4:5) == 'RF' ) then
      write(stderr,*) 'THIS CODE IS NOT SUPPORTED.'
      write(stderr,*) 'CHOSE ONE SCENARIO CODE ',dattyp(1:3),'(26-45-60-85).'
      call die('icbc','Unknown dattyp',1)
    end if
    call init_gn6hnc
  end if

  call newfile(idate)

  do nnn = 1 , nsteps

    if (.not. lsamemonth(idate, iodate) ) then
      call newfile(monfirst(idate))
    end if

    if ( dattyp(1:4) == 'NNRP' .or. dattyp(1:3) == 'CFS' ) then
      call get_ncep(idate)
    else if ( dattyp == 'ECMWF' ) then
      call get_ecwcp(idate)
    else if ( dattyp == 'ERA40' ) then
      call get_era40(idate)
    else if ( dattyp(1:4) == 'ERA5' ) then
      call get_era5(idate)
    else if ( dattyp == 'ERAIN' .or. dattyp(1:3) == 'EIN' .or. &
              dattyp == 'EIXXX' ) then
      call get_ein(idate)
    else if ( dattyp == 'ERAHI' ) then
      call get_ehi(idate)
    else if ( dattyp(1:3) == 'ECE' ) then
      call get_ecens(idate)
    else if ( dattyp(1:2) == 'EH' ) then
      call get_eh5om(idate)
    else if ( dattyp == 'FVGCM' ) then
      call get_fvgcm(idate)
    else if ( dattyp == 'FNEST' ) then
      call get_nest(idate)
    else
      call get_gn6hnc(idate)
    end if

    call writef(idate)

    iodate = idate
    idate = idate + tbdy

  end do

  call close_output
  call closesst

  if ( dattyp(1:4) == 'ERA5' ) then
    call conclude_era5
  else if ( dattyp == 'ERAIN' .or. dattyp(1:3) == 'EIN' .or. &
            dattyp == 'EIXXX' ) then
    call conclude_ein
  else if ( dattyp(1:3) == 'ECE' ) then
    call conclude_ecens
  else if ( dattyp == 'ECMWF' ) then
    call conclude_ecwcp
  else if ( dattyp(1:2) == 'EH' ) then
    call conclude_eh5om
  else if ( dattyp == 'ERA40' ) then
    call conclude_era40
  else if ( dattyp == 'ERAHI' ) then
    call conclude_ehi
  else if ( dattyp == 'FVGCM' ) then
    call conclude_fvgcm
  else if ( dattyp(1:4) == 'NNRP' .or. dattyp(1:3) == 'CFS' ) then
    call conclude_ncep
  else if ( dattyp == 'FNEST' ) then
    call conclude_nest
  else
    call conclude_gn6hnc
  end if

  call pju%destruct( )
  call pjv%destruct( )
  call pjd%destruct( )
  call dispose_output
  call memory_destroy

  call finaltime(0)
  write(stdout,*) 'Successfully completed ICBC'

#ifdef PNETCDF
  call mpi_finalize(ierr)
#endif

end program icbc
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
