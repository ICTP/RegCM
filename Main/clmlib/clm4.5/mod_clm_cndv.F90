module mod_clm_cndv

#if (defined CNDV)
  !
  ! Module containing routines to drive the annual dynamic vegetation
  ! that works with CN, reset related variables,
  ! and initialize/reset time invariant variables
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams
  use mod_mpmessage
  use mod_stdio
  use mod_date
  use mod_mppparam
  use mod_clm_type
  use mod_clm_nchelper
  use mod_clm_time_manager, only : getdatetime
  use mod_clm_varctl, only : caseid, inst_suffix, nextdate
  use mod_clm_varctl, only : ctitle, finidat, fsurdat, fpftcon
  use mod_clm_cnvegstructupdate, only : CNVegStructUpdate
  use mod_clm_cndvestablishment, only : Establishment
  use mod_clm_cndvlight, only : Light
  use mod_clm_decomp, only : get_proc_bounds, get_proc_global, &
         gcomm_gridcell, gcomm_landunit, gcomm_column, gcomm_pft, ldecomp
  use mod_clm_varpar, only : maxpatch_pft
  use mod_clm_domain, only : ldomain
  use mod_clm_varcon, only : spval
  use mod_clm_varcon , only : secspday

  implicit none

  private

  save

  public :: dv              ! Drives the annual dynamic vegetation that
                            ! works with CN
  public :: histCNDV        ! Output CNDV history file

  contains
  !
  ! Drives the annual dynamic vegetation that works with CN
  !
  subroutine dv(lbg, ubg, lbp, ubp, num_natvegp, filter_natvegp, kyr)
    implicit none
    integer(ik4), intent(in) :: lbg, ubg  ! gridcell bounds
    integer(ik4), intent(in) :: lbp, ubp  ! pft bounds
    ! number of naturally-vegetated pfts in filter
    integer(ik4), intent(inout) :: num_natvegp
    ! filter for naturally-vegetated pfts
    integer(ik4), intent(inout) :: filter_natvegp(ubp-lbp+1)
    integer(ik4), intent(in) :: kyr   ! used in routine climate20 below
    integer(ik4), pointer, contiguous :: mxy(:)  ! pft m index (for laixy(i,j,m),etc.)
    integer(ik4), pointer, contiguous :: pgridcell(:)   ! gridcell of corresponding pft
    ! foliar projective cover on gridcell (fraction)
    real(rk8), pointer, contiguous :: fpcgrid(:)
    ! accumulated growing degree days above 5
    real(rk8), pointer, contiguous :: agdd(:)
    ! annual min of t_mo (Kelvin)
    real(rk8), pointer, contiguous :: t_mo_min(:)
    real(rk8), pointer, contiguous :: temp_tmomin(:)
    real(rk8), pointer, contiguous :: temp_agdd(:)
    real(rk8), pointer, contiguous :: temp_count(:)
    real(rk8), pointer, contiguous :: tmomin20(:) ! 20-yr running mean of tmomin
    real(rk8), pointer, contiguous :: agdd20(:)   ! 20-yr running mean of agdd
    integer(ik4) :: g, p             ! indices
    integer(ik4) :: ier

    ! Assign local pointers to derived type members (gridcell-level)

    agdd20    => clm3%g%gdgvs%agdd20
    tmomin20  => clm3%g%gdgvs%tmomin20

    ! Assign local pointers to derived type members (pft-level)

    mxy       => clm3%g%l%c%p%mxy
    pgridcell => clm3%g%l%c%p%gridcell
    fpcgrid   => clm3%g%l%c%p%pdgvs%fpcgrid
    t_mo_min  => clm3%g%l%c%p%pdgvs%t_mo_min
    agdd      => clm3%g%l%c%p%pdgvs%agdd

    ! *************************************************************************
    ! S. Levis version of LPJ's routine climate20: 'Returns' tmomin20 & agdd20
    ! for use in routine bioclim, which I have placed in routine Establishment
    ! Instead of 20-yr running mean of coldest monthly temperature,
    ! use 20-yr running mean of minimum 10-day running mean
    ! *************************************************************************

    allocate(temp_tmomin(lbg:ubg), temp_agdd(lbg:ubg), &
             temp_count(lbg:ubg), stat=ier)
    if (ier /= 0) call fatal(__FILE__,__LINE__,&
      'DV: allocation error for temp_tmomin,temp_agdd,temp_count')

    temp_tmomin(:) = 0.0_rk8
    temp_agdd(:) = 0.0_rk8
    temp_count(:) = 0.0_rk8
    do p = lbp, ubp
       g = pgridcell(p)
       temp_tmomin(g) = temp_tmomin(g) + t_mo_min(p)
       temp_agdd(g)   = temp_agdd(g)   + agdd(p)
       temp_count(g)  = temp_count(g) + 1.0_rk8
    end do

    if ( kyr == 1 ) then
      do g = lbg, ubg
        if ( temp_count(g) > 0.0_rk8 ) then
          tmomin20(g) = temp_tmomin(g)/temp_count(g)
          agdd20(g)   = temp_agdd(g)/temp_count(g)
        end if
      end do
    else
      do g = lbg, ubg
        if ( temp_count(g) > 0.0_rk8 ) then
          tmomin20(g) = (19._rk8 * tmomin20(g) + &
            temp_tmomin(g)/temp_count(g)) / 20._rk8
          agdd20(g)   = (19._rk8 * agdd20(g) +   &
            temp_agdd(g)/temp_count(g)) / 20._rk8
        end if
      end do
    end if

    deallocate(temp_tmomin,temp_agdd,temp_count)

    ! Rebuild filter of present natually-vegetated pfts after Kill()

    call BuildNatVegFilter(lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns fpcgrid and nind

    call Light(lbg, ubg, lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns updated fpcgrid, nind, crownarea, and present. Due to updated
    ! present, we do not use the natveg filter in this subroutine.

    call Establishment(lbg, ubg, lbp, ubp)

    ! Reset dgvm variables needed in next yr (too few to keep subr. dvreset)

    do p = lbp,ubp
       clm3%g%l%c%p%pcs%leafcmax(p) = 0._rk8
       clm3%g%l%c%p%pdgvs%t_mo_min(p) = 1.0e+36_rk8
    end do
  end subroutine dv
  !
  ! Create CNDV history file
  !
  subroutine histCNDV()
    implicit none
    ! true=>landunit is not vegetated (landunit-level)
    logical, pointer, contiguous :: ifspecial(:)
    ! gridcell index of corresponding pft (pft-level)
    integer(ik4), pointer, contiguous :: pgridcell(:)
    ! landunit index of corresponding pft (pft-level)
    integer(ik4), pointer, contiguous :: plandunit(:)
    ! pft m index (for laixy(i,j,m),etc.)
    integer(ik4), pointer, contiguous :: mxy(:)
    ! foliar projective cover on gridcell (fraction)
    real(rk8), pointer, contiguous :: fpcgrid(:)
    ! number of individuals (#/m**2)
    real(rk8), pointer, contiguous :: nind(:)
    character(len=256) :: dgvm_fn   ! dgvm history filename
    type(clm_filetype) :: ncid      ! netcdf file id
    integer(ik4) :: p,l,c         ! indices
    integer(ik4) :: begp, endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc, endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl, endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg, endg ! per-proc gridcell ending gridcell indices
    integer(ik4) :: ihost, iktau
    integer(ik4) :: numg, numl, numc, nump
    integer(ik4) :: ier                     ! error status
    integer(ik4) :: mdcur, mscur, mcdate    ! outputs from curr_time
    integer(ik4) :: yr,mon,day,mcsec        ! outputs from curr_date
    integer(ik4) :: hours,minutes,secs      ! hours,minutes,seconds of hh:mm:ss
    integer(ik4) :: nbsec                   ! seconds components of a date
    real(rk8):: time                   ! current time
    real(rk8), pointer, contiguous, dimension(:) :: rparr     ! temporary
    integer(ik4), pointer, contiguous, dimension(:) :: iparr  ! temporary
    character(len=256) :: str          ! temporary string
    character(len=  8) :: curdate      ! current date
    character(len=  8) :: curtime      ! current time
    character(len= 10) :: basedate     ! base date (yyyymmdd)
    character(len=  8) :: basesec      ! base seconds
    integer(ik4) :: hostnm
    character (len=32) :: hostname='?'
    character (len=32) :: user='?'

    ! Assign local pointers to derived type members (gridcell-level)

    ! NONE

    ! Assign local pointers to derived type members (landunit-level)

    ifspecial  => clm3%g%l%ifspecial

    ! Assign local pointers to derived subtypes components (pft-level)

    mxy       => clm3%g%l%c%p%mxy
    pgridcell => clm3%g%l%c%p%gridcell
    plandunit => clm3%g%l%c%p%landunit
    fpcgrid   => clm3%g%l%c%p%pdgvs%fpcgrid
    nind      => clm3%g%l%c%p%pdgvs%nind

    ! Determine subgrid bounds for this processor and allocate dynamic memory

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg,numl,numc,nump)

    allocate(rparr(begp:endp), stat=ier)
    if (ier /= 0) call fatal(__FILE__,__LINE__,&
      'histCNDV: allocation error for real rparr')
    allocate(iparr(begp:endp), stat=ier)
    if (ier /= 0) call fatal(__FILE__,__LINE__,&
      'histCNDV: allocation error for integer arrays')

    ! -----------------------------------------------------------------------
    ! Create new netCDF file. File will be in define mode
    ! -----------------------------------------------------------------------

    dgvm_fn = set_dgvm_filename()

    call clm_createfile(trim(dgvm_fn),ncid)

    ! -----------------------------------------------------------------------
    ! Create global attributes.
    ! -----------------------------------------------------------------------

    str = 'CF1.0'
    call clm_addatt(ncid, 'conventions', trim(str))

    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call clm_addatt(ncid, 'history', trim(str))

#ifdef IBM
    hostname='ibm platform '
    user= 'Unknown'
#else
    ihost = hostnm(hostname)
    call getlog(user)
#endif

    call clm_addatt(ncid, 'logname', trim(user))
    call clm_addatt(ncid, 'host', trim(hostname))

    str = 'Community Land Model: CLM3'
    call clm_addatt(ncid, 'source',  trim(str))

    str = '$Name$'
    call clm_addatt(ncid, 'version', trim(str))

    str = '$Id$'
    call clm_addatt(ncid, 'revision_id',  trim(str))

    str = ctitle
    call clm_addatt(ncid, 'case_title', trim(str))

    str = caseid
    call clm_addatt(ncid, 'case_id',  trim(str))

    call clm_addatt(ncid, 'Surface_dataset',  trim(fsurdat))

    str = 'arbitrary initialization'
    if (finidat /= ' ') str = finidat
    call clm_addatt(ncid, 'Initial_conditions_dataset',  trim(str))

    call clm_addatt(ncid,'PFT_physiological_constants_dataset',trim(fpftcon))

    ! -----------------------------------------------------------------------
    ! Define dimensions.
    ! -----------------------------------------------------------------------

    call clm_adddim(ncid,'gridcell',numg)
    call clm_adddim(ncid,'landunit',numl)
    call clm_adddim(ncid,'column',numc)
    call clm_adddim(ncid,'pft',nump)
    call clm_adddim(ncid,'jx',njoutsg)
    call clm_adddim(ncid,'iy',nioutsg)
    call clm_adddim(ncid,'time',clmvar_unlim)
    call clm_adddim(ncid,'string_length',80)

    ! -----------------------------------------------------------------------
    ! Define variables
    ! -----------------------------------------------------------------------

    call clm_addvar(clmvar_integer,ncid,'numpft')

    call curr_time(nextdate,mdcur, mscur)
    call ref_date(yr,mon,day,nbsec)
    hours   = nbsec / 3600
    minutes = (nbsec - hours*3600) / 60
    secs    = (nbsec - hours*3600 - minutes*60)
    write(basedate,'(i4.4,"-",i2.2,"-",i2.2)') yr,mon,day
    write(basesec ,'(i2.2,":",i2.2,":",i2.2)') hours, minutes, secs
    str = 'days since ' // basedate // " " // basesec
    time = mdcur + mscur/secspday

    call clm_addvar(clmvar_double,ncid,'time',['time'],'time',trim(str))

    ! Define surface grid (coordinate variables, latitude, longitude,
    ! surface type).

    call clm_addvar(clmvar_double,ncid,'longxy',['gridcell'], &
          long_name='longitude', units='degrees_east')
    call clm_addvar(clmvar_double,ncid,'latixy',['gridcell'], &
          long_name='latitude', units='degrees_north')
    call clm_addvar(clmvar_integer,ncid,'landmask',['gridcell'], &
          long_name='land/ocean mask (0.=ocean and 1.=land)', &
          missing_value=1,fill_value=1)
    call clm_addvar(clmvar_integer,ncid,'pftxgdc',['gridcell'], &
            long_name='Number of pfts per gridcell', &
            missing_value=1,fill_value=1)
    call clm_addvar(clmvar_logical,ncid,'regcm_mask',['jx','iy'], &
          long_name='regcm land mask', units='1')

    ! Define time information

    call clm_addvar(clmvar_integer,ncid,'mcdate',['time'], &
         long_name='current date (YYYYMMDD)')
    call clm_addvar(clmvar_integer,ncid,'mcsec',['time'], &
         long_name='current seconds of current date', units='s')
    call clm_addvar(clmvar_integer,ncid,'mdcur',['time'], &
         long_name='current day (from base day)')
    call clm_addvar(clmvar_integer,ncid,'mscur',['time'], &
         long_name='current seconds of current day', units='s')
    call clm_addvar(clmvar_integer,ncid,'nstep',['time'], &
         long_name='time step', units='s')

    call clm_addvar(clmvar_double,ncid,'pfts1d_wtxy',['pft'], &
            long_name='pft weight relative to corresponding gridcell', &
            missing_value=1,fill_value=1)
    call clm_addvar(clmvar_integer,ncid,'pfts1d_gridcell', &
            ['pft'],long_name='pft gridcell index', &
            missing_value=1,fill_value=1)
    call clm_addvar(clmvar_integer,ncid,'pfts1d_itypveg', &
            ['pft'],long_name='pft vegetation type', &
            missing_value=1,fill_value=1)
    call clm_addvar(clmvar_integer,ncid,'pfts1d_ityplun',['pft'], &
         long_name='pft landunit type (vegetated,urban,lake,wetland,glacier)', &
         missing_value=1,fill_value=1)

    ! Define time dependent variables

     call clm_addvar(clmvar_double,ncid,'FPCGRID', ['pft ','time'], &
          long_name='plant functional type cover', &
          units='% of vegetated area', &
          missing_value=1,fill_value=1)
     call clm_addvar(clmvar_double,ncid,'NIND', ['pft ','time'], &
          long_name='number of individuals', &
          units='individuals/m2 vegetated land', &
          missing_value=1,fill_value=1)

    call clm_enddef(ncid)

    ! -----------------------------------------------------------------------
    ! Write variables
    ! -----------------------------------------------------------------------

    call clm_writevar(ncid,'numpft',maxpatch_pft)

    ! Write surface grid (coordinate variables, latitude,
    ! longitude, surface type).

    call clm_writevar(ncid,'longxy',ldomain%lonc,gcomm_gridcell)
    call clm_writevar(ncid,'latixy',ldomain%latc,gcomm_gridcell)
    call clm_writevar(ncid,'landmask',ldomain%mask,gcomm_gridcell)
    call clm_writevar(ncid,'pftxgdc',ldecomp%pftxgdc,gcomm_gridcell)
    call clm_writevar(ncid,'regcm_mask',lndcomm%global_out_sgmask)
    call clm_writevar(ncid,'pfts1d_wtxy',clm3%g%l%c%p%wtgcell,gcomm_pft)
    do p = begp, endp
      iparr(p) = clm3%g%l%c%p%gridcell(p)
    end do
    call clm_writevar(ncid,'pfts1d_gridcell',iparr,gcomm_pft)
    call clm_writevar(ncid,'pfts1d_itypveg',clm3%g%l%c%p%itype,gcomm_pft)
    do p = begp, endp
      iparr(p) = clm3%g%l%itype(clm3%g%l%c%p%landunit(p))
    end do
    call clm_writevar(ncid,'pfts1d_ityplun',iparr,gcomm_pft)

    ! Write current date, current seconds, current day, current nstep

    call curr_date(nextdate,yr,mon,day,mcsec)
    mcdate = yr*10000 + mon*100 + day

    iktau = int(syncro_srf%lcount,ik4)
    call clm_writevar(ncid,'mcdate',mcdate,nt=1)
    call clm_writevar(ncid,'mcsec',mcsec,nt=1)
    call clm_writevar(ncid,'mdcur',mdcur,nt=1)
    call clm_writevar(ncid,'mscur',mscur,nt=1)
    call clm_writevar(ncid,'nstep',iktau,nt=1)
    call clm_writevar(ncid,'time',time,nt=1)

    ! Write time dependent variables to CNDV history file

    ! The if .not. ifspecial statment below guarantees that the m index will
    ! always lie between 1 and maxpatch_pft

    rparr(:) = 0._rk8
    do p = begp, endp
       l = plandunit(p)
       if (.not. ifspecial(l)) rparr(p) = fpcgrid(p)*100._rk8
    end do
    call clm_writevar(ncid,'FPCGRID',rparr,gcomm_pft,nt=1)

    rparr(:) = 0._rk8
    do p = begp, endp
       l = plandunit(p)
       if (.not. ifspecial(l)) rparr(p) = nind(p)
    end do
    call clm_writevar(ncid,'NIND',rparr,gcomm_pft,nt=1)

    ! Deallocate dynamic memory

    deallocate(rparr)
    deallocate(iparr)

    !------------------------------------------------------------------
    ! Close and archive netcdf CNDV history file
    !------------------------------------------------------------------

    call clm_closefile(ncid)

    if (myid == italk) then
       write(stdout,*) 'Written CNDV history dataset'
    end if

  end subroutine histCNDV
  !
  ! Determine initial dataset filenames
  !
  character(len=256) function set_dgvm_filename ()
    implicit none
    character(len=4) :: cdate  !date char string
    integer(ik4) :: day        !day (1 -> 31)
    integer(ik4) :: mon        !month (1 -> 12)
    integer(ik4) :: yr         !year (0 -> ...)
    integer(ik4) :: sec        !seconds into current day
    call curr_date(nextdate,yr,mon,day,sec)
    write(cdate,'(i4.4)') yr+1
    set_dgvm_filename = trim(dirout)//pthsep//trim(caseid)//  &
                              ".clm."//trim(inst_suffix)// &
                              ".hv."// trim(cdate) //".nc"
  end function set_dgvm_filename
  !
  ! Reconstruct a filter of naturally-vegetated PFTs for use in DGVM
  !
  subroutine BuildNatVegFilter(lbp, ubp, num_natvegp, filter_natvegp)
    implicit none
    ! pft bounds
    integer(ik4), intent(in)  :: lbp, ubp
    ! number of pfts in naturally-vegetated filter
    integer(ik4), intent(out) :: num_natvegp
    ! pft filter for naturally-vegetated points
    integer(ik4), intent(out) :: filter_natvegp(ubp-lbp+1)
    logical, pointer, contiguous :: present(:)     ! whether this pft present in patch
    integer(ik4) :: p

    ! Assign local pointers to derived type members (pft-level)
    present   => clm3%g%l%c%p%pdgvs%present

    num_natvegp = 0
    do p = lbp,ubp
      if (present(p)) then
        num_natvegp = num_natvegp + 1
        filter_natvegp(num_natvegp) = p
      end if
    end do
  end subroutine BuildNatVegFilter

#endif

end module mod_clm_cndv
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
