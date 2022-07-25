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

module mod_clm_params

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_stdio
  use mod_dynparam
  use mod_date
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_domain
  use mod_service
  use mod_sun
  use mod_atm_stub
  use mod_zita
#ifdef CLM45
  use mod_clm_regcm
#endif
  use mod_ipcc_scenario , only : set_scenario
  use mod_ncio
  use mod_timer

  implicit none

  private

  real(rkx) , parameter :: mindt = 1.0_rkx

  public :: param

  contains
  !
  ! This subroutine defines the various model parameters.
  !
  subroutine param
    implicit none
    integer(ik4) :: iretval
    integer(ik4) :: i , j , k
    integer(ik8) :: mdate0 , mdate1 , mdate2
    integer(ik4) :: hspan , ipunit
    integer(ik4) :: n , len_path
    character(len=32) :: appdat
    type(rcm_time_interval) :: bdif
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'param'
    integer(ik4) , save :: idindx = 0
#endif
    !
    ! namelist:
    !
    namelist /restartparam/ ifrest , mdate0 , mdate1 , mdate2

    namelist /timeparam/ dtrad , dtsrf , dtcum , dtche , dtabem , dt

    namelist /tweakparam/ itweak_temperature , itweak_solar_irradiance , &
            itweak_sst , itweak_greenhouse_gases , temperature_tweak ,   &
            sst_tweak , solar_tweak , gas_tweak_factors

    namelist /clmsaparam/ dirout , sfbcread

#ifdef DEBUG
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! default values for all the options:
    !     (can be overwritten by namelist input).
    !
    ! restartparam ;
    !
    ifrest = .false.     ! Restart?:  t=true; f=false
    idate0 = 1900010100  ! Start date of the initial simulation
    idate1 = 1900010100  ! Start date of this simulation
    idate2 = 1900010100  ! End Date this simulation
    !
    ! timeparam ;
    !
    dt = 100.0_rkx  ! time step in seconds
    dtrad = 0.0_rkx ! time interval in min solar rad caluclated
    dtsrf = 0.0_rkx ! time interval at which bats is called (secs)
    dtcum = 0.0_rkx ! time interval at which cumulus is called (secs)
    dtabem = 0.0_rkx ! time interval absorption-emission calculated (hours)
    dtche = 900.0_rkx ! time interval at which bats is called (secs)
    dirout = './output'
    sfbcread = 1
    lsync = .true.
    do_parallel_netcdf_in = .false.
    scenario = 'RCP4.5'
    ghg_year_const = 1950
    ifixsolar = 0
    isolconst = 1
    fixedsolarval = 1367.0_rkx
    !
    ! tweakparam ;
    !
    itweak = 0
    itweak_sst = 0
    itweak_temperature = 0
    itweak_solar_irradiance = 0
    itweak_greenhouse_gases = 0
    sst_tweak = 0.0_rkx
    temperature_tweak = 0.0_rkx
    solar_tweak = 0.0_rkx
    gas_tweak_factors(:) = 1.0_rkx

#ifdef CLM
    if ( myid == italk ) then
      if (nsg /= 1 ) then
        write (stderr,*) 'Running SUBGRID with CLM: not implemented'
        write (stderr,*) 'Please set nsg to 1 in regcm.in'
        call fatal(__FILE__,__LINE__, &
                   'CLM & SUBGRID TOGETHER')
      end if
    end if
#endif

    if ( myid == iocpu ) then
      open(newunit=ipunit, file=namelistfile, status='old', &
                   action='read', iostat=iretval)
      if ( iretval /= 0 ) then
        write(stderr,*) 'Error opening input namelist file ',trim(namelistfile)
        call fatal(__FILE__,__LINE__, &
                   'INPUT NAMELIST OPEN ERROR')
#ifdef DEBUG
      else
        write(stdout,*) 'Open ',trim(namelistfile),' OK'
#endif
      end if

      write(stdout,*) 'Reading model namelist file'

      rewind(ipunit)
      read (ipunit, nml=restartparam, iostat=iretval, err=100)
      if ( iretval /= 0 ) then
        write(stderr,*) 'Error reading restartparam namelist'
        call fatal(__FILE__,__LINE__, &
                   'INPUT NAMELIST READ ERROR')
#ifdef DEBUG
      else
        write(stdout,*) 'Read restartparam OK'
#endif
      end if

      rewind(ipunit)
      read (ipunit, nml=clmsaparam, iostat=iretval, err=200)
      if ( iretval /= 0 ) then
        write(stderr,*) 'Error reading clmsaparam namelist'
        call fatal(__FILE__,__LINE__, &
                   'INPUT NAMELIST READ ERROR')
#ifdef DEBUG
      else
        write(stdout,*) 'Read clmsaparam OK'
#endif
      end if

      idate0 = i8wcal(mdate0,ical)
      idate1 = i8wcal(mdate1,ical)
      idate2 = i8wcal(mdate2,ical)
      bdif = idate2 - idate1
      hspan = nint(tohours(bdif))
      if ( mod(hspan,24) /= 0 ) then
        call fatal(__FILE__,__LINE__,  &
                   'Runtime increments must be modulus 24 hours')
      end if

      rewind(ipunit)
      read (ipunit, nml=timeparam, iostat=iretval, err=101)
      if ( iretval /= 0 ) then
        write(stderr,*) 'Error reading timeparam namelist'
        call fatal(__FILE__,__LINE__, &
                   'INPUT NAMELIST READ ERROR')
#ifdef DEBUG
      else
        write(stdout,*) 'Read timeparam OK'
#endif
      end if

      len_path = len(trim(dirout))
      if ( dirout(len_path:len_path) /= '/' ) dirout = trim(dirout)//'/'

      if ( itweak == 1 ) then
        rewind(ipunit)
        read (ipunit , tweakparam, iostat=iretval, err=121)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Tweak parameters absent.'
          write(stdout,*) 'Disable tweaking.'
          itweak = 0
#ifdef DEBUG
        else
          write(stdout,*) 'Read tweakparam OK'
#endif
        end if
        if ( itweak_sst == 0 .and.              &
             itweak_temperature == 0 .and.      &
             itweak_solar_irradiance == 0 .and. &
             itweak_greenhouse_gases == 0 ) then
          write(stdout,*) 'Tweak parameters not enabled.'
          write(stdout,*) 'Disable tweaking.'
          itweak = 0
        end if
      end if

      close(ipunit)

      if ( dt < mindt ) then
        write(stderr,*) 'Minimum dt allowed is ',mindt,' seconds.'
        call fatal(__FILE__,__LINE__, &
                   'DT TOO SMALL')
      end if

      dt = check_against_outparams(dt,mindt)

      if ( dtsrf <= 0.0_rkx ) dtsrf = 600.0_rkx
      if ( dtcum <= 0.0_rkx ) dtcum = 300.0_rkx
      if ( dtche <= 0.0_rkx ) dtche = 900.0_rkx
      if ( dtrad <= 0.0_rkx ) then
        dtrad = 1800.0_rkx
      else
        dtrad = dtrad * 60.0_rkx
      end if
      if ( dtabem <= 0.0_rkx ) then
        dtabem = 64800.0_rkx
      else
        dtabem = dtabem * 3600.0_rkx
      end if

      if ( dtcum < dt ) dtcum = dt
      if ( dtsrf < dt ) dtsrf = dt
      if ( dtche < dt ) dtche = dt
      if ( dtrad < dt ) dtrad = dt
      if ( dtabem < dt ) dtabem = dt

      dtsrf = int(dtsrf / dt) * dt
      if ( ifshf ) then
        do while ( mod(3600.0_rkx,dtsrf) > d_zero )
          dtsrf = dtsrf - dt
        end do
      else
        do while ( mod(86400.0_rkx,dtsrf) > d_zero )
          dtsrf = dtsrf - dt
        end do
      end if
      dtsrf = int(dtsrf / dt) * dt

      dtcum = int(dtcum / dt) * dt
      do while ( mod(86400.0_rkx,dtcum) > d_zero )
        dtcum = dtcum - dt
      end do
      dtcum = int(dtcum / dt) * dt

      dtrad = int(dtrad / dt) * dt
      do while ( mod(86400.0_rkx,dtrad) > d_zero )
        dtrad = dtrad - dt
      end do
      dtrad = int(dtrad / dt) * dt

      dtche = int(dtche / dt) * dt
      dtabem = int(dtabem / dtrad) * dtrad

      if ( iseaice == 1 ) then
        select case (ssttyp)
          case ('EIN15','EIN75','EIXXX','ERA5 ')
            icetriggert = 271.465_rkx
          case default
            icetriggert = 271.355_rkx
        end select
      end if
    end if
    !
    ! communicate to all processors
    !
    call bcast(lsmoist)
    call bcast(ifrest)
    call bcast(hspan)
    call bcast(idate0)
    call bcast(idate1)
    call bcast(idate2)
    call bcast(globidate1)
    call bcast(globidate2)

    call bcast(dt)
    call bcast(dtrad)
    call bcast(dtsrf)
    call bcast(dtcum)
    call bcast(dtche)
    call bcast(dtabem)
    call bcast(itweak)

    call bcast(sfbcread)

    if ( itweak == 1 ) then
      call bcast(itweak_sst)
      call bcast(itweak_temperature)
      call bcast(itweak_solar_irradiance)
      call bcast(itweak_greenhouse_gases)
      call bcast(sst_tweak)
      call bcast(temperature_tweak)
      call bcast(solar_tweak)
      call bcast(gas_tweak_factors)
    end if

    rcmtimer => rcm_timer(idate0,idate1,idate2,dtsrf)

    !
    ! ALLOCATE NEEDED SPACE
    !

    call allocate_mod_runparams
    call allocate_mod_atm_interface

    if ( ifrest ) then
      doing_restart = .true.
    end if
    !
    ! Calculate the time step in minutes.
    !
    dtsec = dt
    rdt   = d_one/dt
    dtbdys = real(ibdyfrq,rkx)*secph
    !
    ! Reset the options/calculate variables using namelist info
    !
    bdydate1 = idate1

    alarm_hour => rcm_alarm(rcmtimer,3600.0_rkx)
    alarm_day => rcm_alarm(rcmtimer,86400.0_rkx)
    alarm_in_bdy => rcm_alarm(rcmtimer,dtbdys)

    syncro_rep => rcm_syncro(rcmtimer,3.0_rkx*3600.0_rkx)
    syncro_srf => rcm_syncro(rcmtimer,dtsrf)
    syncro_cum => rcm_syncro(rcmtimer,dtcum)
    syncro_rad => rcm_syncro(rcmtimer,dtrad)
    syncro_emi => rcm_syncro(rcmtimer,dtabem)
    syncro_che => rcm_syncro(rcmtimer,dtche)
    if ( debug_level > 0 ) then
      syncro_dbg => rcm_syncro(rcmtimer,secph*dbgfrq)
    end if

    rnsrf_for_day = syncro_srf/alarm_day

    dtsq = dt*dt
    dtcb = dt*dt*dt

    intbdy = rcm_time_interval(ibdyfrq,uhrs)
    intsom = rcm_time_interval(1,umnt)

    if ( myid == italk ) then
      appdat = tochar(idate0)
      write(stdout,*) 'Initial date of the global simulation: ', appdat
      appdat = tochar(idate1)
      write(stdout,*) 'Initial date of this run             : ', appdat
      appdat = tochar(idate2)
      write(stdout,*) 'Final date of this run               : ', appdat
      write(stdout,*) 'Total simulation lenght              : ', hspan, ' hours'
      write(stdout,'(a,f11.6)') ' Timestep in seconds = ', dtsec
    end if

    call bcast(dirter,256)
    call bcast(dirglob,256)
    call bcast(dirout,256)
    call bcast(domname,64)
    call read_domain_info(mddom%ht,mddom%lndcat,mddom%lndtex,            &
                          mddom%mask,mddom%area,                         &
                          mddom%xlat,mddom%xlon,mddom%dlat,mddom%dlon,   &
                          mddom%ulat,mddom%ulon,mddom%vlat,mddom%vlon,   &
                          mddom%msfx,mddom%msfd,mddom%msfu,mddom%msfv,   &
                          mddom%coriol,mddom%snowam,mddom%smoist,        &
                          mddom%rmoist,mddom%dhlake,base_state_ts0)
    call bcast(ds)
    call bcast(ptop)
    call bcast(xcone)

    dx = ds * d_1000
    dx2 = d_two*dx
    dx4 = d_four*dx
    dx8 = 8.0_rkx*dx
    dx16 = 16.0_rkx*dx
    dxsq = dx*dx
    rdx = 1.0_rkx/dx
    rdxsq = 1.0_rkx/dxsq

    if ( myid == italk ) then
      write(stdout,*) 'Setting IPCC scenario to ', scenario
      if ( scenario == 'CONST' ) then
        write(stdout,*) 'Selected value at year ', ghg_year_const
      end if
    end if

    call set_scenario(scenario,rcmtimer%year,rcmtimer%month)

    if ( myid == italk ) then
      if ( ifrest .and. idate0 == idate1 ) then
        write(stderr,*) 'Error in parameter set.'
        write(stderr,*) 'Cannot set idate0 == idate1 on restart run'
        write(stderr,*) 'Correct idate0.'
        call fatal(__FILE__,__LINE__, &
                   'IDATE0==IDATE1 ON RESTART')
      else if ( .not. ifrest .and. idate0 /= idate1 ) then
        write(stderr,*) 'Error in parameter set.'
        write(stderr,*) 'Cannot set idate0 /= idate1 on non restart run'
        write(stderr,*) 'Correct idate1.'
        call fatal(__FILE__,__LINE__, &
                   'IDATE0/=IDATE1 ON NON RESTART')
      end if
    end if

    if ( nsg > 1 ) then
      call read_subdomain_info(mdsub%ht,mdsub%lndcat,mdsub%lndtex,mdsub%mask, &
               mdsub%area,mdsub%xlat,mdsub%xlon,mdsub%dhlake)
      mdsub%ht = mdsub%ht*egrav
    else
      do i = ici1 , ici2
        do j = jci1 , jci2
          mdsub%ht(1,j,i) = mddom%ht(j,i)*egrav
          mdsub%lndcat(1,j,i) = mddom%lndcat(j,i)
          mdsub%lndtex(1,j,i) = mddom%lndtex(j,i)
          mdsub%xlat(1,j,i) = mddom%xlat(j,i)
          mdsub%xlon(1,j,i) = mddom%xlon(j,i)
          mdsub%mask(1,j,i) = mddom%mask(j,i)
          mdsub%area(1,j,i) = mddom%area(j,i)
        end do
      end do
    end if
    !
    !------invert mapscale factors and convert hgt to geopotential
    !
    do i = ide1 , ide2
      do j = jde1 , jde2
        mddom%ht(j,i)   = mddom%ht(j,i)*egrav
      end do
    end do
    if ( idynamic < 3 ) then
      do k = 1 , kz
        dsigma(k) = (sigma(k+1) - sigma(k))
        hsigma(k) = (sigma(k+1) + sigma(k))*d_half
      end do
    else
      call compute_moloch_static
    end if
    !
    !-----compute land/water mask
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( mddom%mask(j,i) > 0.1_rkx ) then
          mddom%ldmsk(j,i) = 1
        else
          mddom%ldmsk(j,i) = 0
        end if
      end do
    end do
    !
    !-----compute land/water mask on subgrid space
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( mdsub%mask(n,j,i) > 0.1_rkx ) then
            mdsub%ldmsk(n,j,i) = 1
          else
            mdsub%ldmsk(n,j,i) = 0
          end if
        end do
      end do
    end do

    if ( myid == italk ) then
      write(stdout,*) 'Domain grid parameters:'
      write(stdout,'(a,a)') '  Map Projection        : ',iproj
      write(stdout,'(a,i4,a,i4,a,i3)') &
        '  Dot Grid Full Extent  : ',jx,'x',iy,'x',kz
      write(stdout,'(a,f11.6,a)') '  Model Top Pressure    : ',ptop,' cb'
      write(stdout,'(a,f11.6,a)') '  Model Grid Spacing    : ',ds,' km'
      write(stdout,'(a,f11.6,a)') '  Proj Center Latitude  : ',clat,' deg'
      write(stdout,'(a,f11.6,a)') '  Proj Center longitude : ',clon,' deg'
      if ( iproj == 'ROTMER' .or. iproj == 'ROTLLR' ) then
        write(stdout,'(a,f11.6,a)') '  Pole Latitude         : ',plat,' deg'
        write(stdout,'(a,f11.6,a)') '  Pole longitude        : ',plon,' deg'
      else if ( iproj == 'LAMCON' ) then
        write(stdout,'(a,f11.6,a)') '  True Latitude 1       : ',truelatl,' deg'
        write(stdout,'(a,f11.6,a)') '  True Latitude 2       : ',truelath,' deg'
      end if
    end if

    if (itweak == 1 ) then
      if ( myid == italk ) then
        write(stdout,*) 'TWEAKING OF DATA ENABLED!'
        write(stdout,*) 'THIS RUN IS TO BE CONSIDERED A NON STANDARD SCENARIO!'
        if ( itweak_temperature == 1 ) then
          write(stdout,'(a,f11.6,a)') ' Value added to temperature      : ', &
                  temperature_tweak , ' K'
        end if
        if ( itweak_solar_irradiance == 1 ) then
          write(stdout,'(a,f11.6,a)') ' Value added to solar irradiance : ', &
                  solar_tweak , ' W m-2'
        end if
        if ( itweak_greenhouse_gases == 1 ) then
          write(stdout,'(a,f11.6)') ' CO2 concentration factor        : ', &
                  gas_tweak_factors(1)
          write(stdout,'(a,f11.6)') ' CH4 concentration factor        : ', &
                  gas_tweak_factors(2)
          write(stdout,'(a,f11.6)') ' N2O concentration factor        : ', &
                  gas_tweak_factors(3)
          write(stdout,'(a,f11.6)') ' CFC11 concentration factor      : ', &
                  gas_tweak_factors(4)
          write(stdout,'(a,f11.6)') ' CFC12 concentration factor      : ', &
                  gas_tweak_factors(5)
        end if
      end if
    end if

    do i = ici1 , ici2
      do j = jci1 , jci2
        mddom%iveg(j,i) = nint(mddom%lndcat(j,i))
        mddom%itex(j,i) = nint(mddom%lndtex(j,i))
        do n = 1 , nnsg
          mdsub%iveg(n,j,i) = nint(mdsub%lndcat(n,j,i))
          mdsub%itex(n,j,i) = nint(mdsub%lndtex(n,j,i))
        end do
      end do
    end do

    call allocate_surface_model
    call init_surface_model
    call zenitm(mddom%xlat,mddom%xlon,coszrs)

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
    return

100 call fatal(__FILE__,__LINE__, 'Error reading RESTARTPARAM')
101 call fatal(__FILE__,__LINE__, 'Error reading TIMEPARAM')
121 call fatal(__FILE__,__LINE__, 'Error reading TWEAKPARAM')
200 call fatal(__FILE__,__LINE__, 'Error reading CLMSAPARAM')

    contains

      subroutine init_surface_model
        implicit none
      end subroutine init_surface_model

      recursive integer function gcd_rec(u,v) result(gcd)
        implicit none
        integer , intent(in) :: u , v
        if ( mod(u,v) /= 0 ) then
          gcd = gcd_rec(v,mod(u,v))
        else
          gcd = v
        end if
      end function gcd_rec

      real(rkx) function check_against_outparams(dt,dec) result(newdt)
        implicit none
        real(rkx) , intent(in) :: dt , dec
        newdt = int(dt/dec)*dec
        if ( ifshf ) then
          do
            if ( gcd_rec(int(newdt), int(secph)) < newdt ) then
              newdt = newdt + dec
              cycle
            end if
            exit
          end do
        else
          do
            if ( gcd_rec(int(newdt), int(secpd)) < newdt ) then
              newdt = newdt + dec
              cycle
            end if
            exit
          end do
        end if
      end function check_against_outparams

      subroutine compute_moloch_static
        implicit none
        integer(ik4) :: i , j
        real(rkx) , dimension(kzp1) :: fak , fbk
        call model_zitaf(zita)
        call model_zitah(zitah)
        mo_dzita = zita(kz)
        sigma = sigmazita(zita)
        hsigma = sigmazita(zitah)
        fak = md_ak(zita)
        fbk = md_bk(zita)
        ak = md_ak(zitah)
        bk = md_bk(zitah)
        do k = 1 , kz
          dsigma(k) = (sigma(k+1)-sigma(k))
        end do
        do i = ice1 , ice2
          do j = jce1 , jce2
            zeta(j,i) = ak(kz) + (bk(kz) - d_one) * mddom%ht(j,i)*regrav
            fmzf(j,i) = md_fmz(zita(kzp1),mddom%ht(j,i))
          end do
        end do
      end subroutine compute_moloch_static

  end subroutine param

end module mod_clm_params

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
