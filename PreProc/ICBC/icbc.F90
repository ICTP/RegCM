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

subroutine myabort
  implicit none
  call abort
end subroutine myabort

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
  use mod_erahi
  use mod_ecens
  use mod_fvgcm
  use mod_ncep
  use mod_nest
  use mod_gn6hnc
  use mod_write

  implicit none

  integer(ik4) :: nnn
  type(rcm_time_and_date) :: idate , iodate
  type(rcm_time_interval) :: tdiff , tbdy
  integer(ik4) :: nsteps
  integer(ik4) :: ierr
  integer(ik4) :: ifupr
  character(len=256) :: namelistfile, prgname

  namelist /nonhydroparam/ ifupr , logp_lrate

  call header('icbc')
  !
  ! Read input global namelist
  !
  call get_command_argument(0,value=prgname)
  call get_command_argument(1,value=namelistfile)
  call initparam(namelistfile, ierr)
  if ( idynamic == 2 ) then
    open(ipunit, file=namelistfile, status='old', &
         action='read', iostat=ierr)
    rewind(ipunit)
    read(ipunit, nml=nonhydroparam, iostat=ierr)
    if ( ierr /= 0 ) then
      write(stdout, *) 'Using default non hydrostatic parameters'
      ierr = 0
    end if
    close(ipunit)
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

  if (dattyp == 'CCSMN' .or. dattyp == 'CAM4N' .or. dattyp == 'CCSM3' .or. &
      dattyp(1:3) == 'CA_' .or. dattyp(1:3) == 'GF_' ) then
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
  if (dattyp(1:3) == 'MP_' .or. dattyp(1:2) == 'E5' ) then
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
    call headernc
  else if ( dattyp == 'ECMWF' ) then
    call headerec
  else if ( dattyp == 'ERA40' ) then
    call headerera
  else if ( dattyp == 'ERAIN' .or. dattyp(1:3) == 'EIN' .or. &
            dattyp == 'EIXXX' ) then
    call headerein
  else if ( dattyp(1:3) == 'ECE' ) then
    call headerecens
  else if ( dattyp == 'GFS11' ) then
    call headgn6hnc
  else if ( dattyp == 'ERAHI' ) then
    call headerehi
  else if ( dattyp(1:2) == 'EH' ) then
    call headermpi
  else if ( dattyp == 'FVGCM' ) then
    call headerfv
  else if ( dattyp == 'FNEST' ) then
    call headernest
  else if ( dattyp == 'CAM4N'    .or. dattyp == 'CCSMN'    .or. &
            dattyp == 'CCSM3'    .or. dattyp == 'JRA55'    .or. &
            dattyp(1:3) == 'HA_' .or. dattyp(1:3) == 'CA_' .or. &
            dattyp(1:3) == 'IP_' .or. dattyp(1:3) == 'EC_' .or. &
            dattyp(1:3) == 'GF_' .or. dattyp(1:3) == 'CN_' .or. &
            dattyp(1:3) == 'CS_' .or. dattyp(1:3) == 'MP_' .or. &
            dattyp(1:3) == 'MI_' .or. dattyp(1:2) == 'E5' ) then
    call headgn6hnc
  else
    call die('icbc','Unknown dattyp',1)
  end if

  call newfile(idate)

  do nnn = 1 , nsteps

    if (.not. lsamemonth(idate, iodate) ) then
      call newfile(monfirst(idate))
    end if

    if ( dattyp(1:4) == 'NNRP' .or. dattyp(1:3) == 'CFS' ) then
      call getncep(idate)
    else if ( dattyp == 'ECMWF' ) then
      call getecwcp(idate)
    else if ( dattyp == 'ERA40' ) then
      call getera40(idate)
    else if ( dattyp == 'ERAIN' .or. dattyp(1:3) == 'EIN' .or. &
              dattyp == 'EIXXX' ) then
      call getein(idate)
    else if ( dattyp == 'GFS11' ) then
      call get_gn6hnc(idate)
    else if ( dattyp == 'ERAHI' ) then
      call geterahi(idate)
    else if ( dattyp(1:3) == 'ECE' ) then
      call getecens(idate)
    else if ( dattyp(1:2) == 'EH' ) then
      call geteh5om(idate)
    else if ( dattyp == 'FVGCM' ) then
      call getfvgcm(idate)
    else if ( dattyp == 'FNEST' ) then
      call get_nest(idate)
    else if ( dattyp == 'CAM4N'    .or. dattyp == 'CCSMN'    .or. &
              dattyp == 'CCSM3'    .or. dattyp == 'JRA55'    .or. &
              dattyp(1:3) == 'HA_' .or. dattyp(1:3) == 'CA_' .or. &
              dattyp(1:3) == 'IP_' .or. dattyp(1:3) == 'EC_' .or. &
              dattyp(1:3) == 'GF_' .or. dattyp(1:3) == 'CN_' .or. &
              dattyp(1:3) == 'CS_' .or. dattyp(1:3) == 'MP_' .or. &
              dattyp(1:3) == 'MI_' .or. dattyp(1:2) == 'E5' ) then
      call get_gn6hnc(idate)
    end if
    call writef(idate)

    iodate = idate
    idate = idate + tbdy

  end do

  call close_output
  call closesst

  call dispose_output
  call memory_destroy

  call finaltime(0)
  write(stdout,*) 'Successfully completed ICBC'

end program icbc
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
