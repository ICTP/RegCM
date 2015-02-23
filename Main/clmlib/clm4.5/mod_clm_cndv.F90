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
  use mod_date
  use mod_mppparam
  use mod_clm_type
  use mod_clm_nchelper
  use mod_clm_time_manager , only : getdatetime
  use mod_clm_varctl , only : caseid , inst_suffix
  use mod_clm_varctl , only : ctitle, finidat, fsurdat, fpftcon
  use mod_clm_cnvegstructupdate , only : CNVegStructUpdate
  use mod_clm_cndvestablishment , only : Establishment
  use mod_clm_cndvlight , only : Light
  use mod_clm_decomp , only : get_proc_bounds , get_proc_global
  use mod_clm_varpar , only : maxpatch_pft
  use mod_clm_domain , only : ldomain
  use mod_clm_varcon , only : spval
  use mod_clm_varcon  , only : secspday

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
    integer(ik4), intent(in) :: lbg , ubg  ! gridcell bounds
    integer(ik4), intent(in) :: lbp , ubp  ! pft bounds
    ! number of naturally-vegetated pfts in filter
    integer(ik4), intent(inout) :: num_natvegp
    ! filter for naturally-vegetated pfts
    integer(ik4), intent(inout) :: filter_natvegp(ubp-lbp+1)
    integer(ik4), intent(in) :: kyr   ! used in routine climate20 below
    integer(ik4) , pointer :: mxy(:)  ! pft m index (for laixy(i,j,m),etc.)
    integer(ik4) , pointer :: pgridcell(:)   ! gridcell of corresponding pft
    ! foliar projective cover on gridcell (fraction)
    real(rk8), pointer :: fpcgrid(:)
    ! accumulated growing degree days above 5
    real(rk8), pointer :: agdd(:)
    ! annual min of t_mo (Kelvin)
    real(rk8), pointer :: t_mo_min(:)
    real(rk8), pointer :: tmomin20(:) ! 20-yr running mean of tmomin
    real(rk8), pointer :: agdd20(:)   ! 20-yr running mean of agdd
    integer(ik4)  :: g , p            ! indices

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

    do p = lbp,ubp
       g = pgridcell(p)
       if (kyr == 2) then ! slevis: add ".and. start_type==arb_ic" here?
          tmomin20(g) = t_mo_min(p) ! NO, b/c want to be able to start dgvm
          agdd20(g) = agdd(p)       ! w/ clmi file from non-dgvm simulation
       end if
       tmomin20(g) = (19.D0 * tmomin20(g) + t_mo_min(p)) / 20.D0
       agdd20(g)   = (19.D0 * agdd20(g)   + agdd(p)    ) / 20.D0
    end do

    ! Rebuild filter of present natually-vegetated pfts after Kill()

    call BuildNatVegFilter(lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns fpcgrid and nind

    call Light(lbg, ubg, lbp, ubp, num_natvegp, filter_natvegp)

    ! Returns updated fpcgrid, nind, crownarea, and present. Due to updated
    ! present, we do not use the natveg filter in this subroutine.

    call Establishment(lbg, ubg, lbp, ubp)

    ! Reset dgvm variables needed in next yr (too few to keep subr. dvreset)

    do p = lbp,ubp
       clm3%g%l%c%p%pcs%leafcmax(p) = 0.D0
       clm3%g%l%c%p%pdgvs%t_mo_min(p) = 1.0D+36
    end do
  end subroutine dv
  !
  ! Create CNDV history file
  !
  subroutine histCNDV()
    implicit none
    ! true=>landunit is not vegetated (landunit-level)
    logical , pointer :: ifspecial(:)
    ! gridcell index of corresponding pft (pft-level)
    integer(ik4) , pointer :: pgridcell(:)
    ! landunit index of corresponding pft (pft-level)
    integer(ik4) , pointer :: plandunit(:)
    ! pft m index (for laixy(i,j,m),etc.)
    integer(ik4) , pointer :: mxy(:)
    ! foliar projective cover on gridcell (fraction)
    real(rk8), pointer :: fpcgrid(:)
    ! number of individuals (#/m**2)
    real(rk8), pointer :: nind(:)
    character(len=256) :: dgvm_fn   ! dgvm history filename
    type(clm_filetype) :: ncid      ! netcdf file id
    integer(ik4) :: g,p,l           ! indices
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    integer(ik4) :: ihost , iktau
    integer(ik4) :: ier                     ! error status
    integer(ik4) :: mdcur, mscur, mcdate    ! outputs from curr_time
    integer(ik4) :: yr,mon,day,mcsec        ! outputs from curr_date
    integer(ik4) :: hours,minutes,secs      ! hours,minutes,seconds of hh:mm:ss
    integer(ik4) :: nbsec                   ! seconds components of a date
    real(rk8):: time                   ! current time
    character(len=256) :: str          ! temporary string
    character(len=  8) :: curdate      ! current date
    character(len=  8) :: curtime      ! current time
    character(len= 10) :: basedate     ! base date (yyyymmdd)
    character(len=  8) :: basesec      ! base seconds
    real(rk8), pointer :: rbuf2dg(:)  ! temporary
    character(len=32) :: subname='histCNDV'
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

    allocate(rbuf2dg(begp,endp), stat=ier)
    if (ier /= 0) call fatal(__FILE__,__LINE__,&
      'histCNDV: allocation error for rbuf2dg')

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

    call clm_adddim(ncid,'gridcell',ldomain%ns)
    call clm_adddim(ncid,'pft',maxpatch_pft)
    call clm_adddim(ncid,'time',clmvar_unlim)
    call clm_adddim(ncid,'string_length',80)

    ! -----------------------------------------------------------------------
    ! Define variables
    ! -----------------------------------------------------------------------

    call curr_time(idatex,mdcur, mscur)
    call ref_date(yr,mon,day,nbsec)
    hours   = nbsec / 3600
    minutes = (nbsec - hours*3600) / 60
    secs    = (nbsec - hours*3600 - minutes*60)
    write(basedate,'(i4.4,'-',i2.2,'-',i2.2)') yr,mon,day
    write(basesec ,'(i2.2,':',i2.2,':',i2.2)') hours, minutes, secs
    str = 'days since ' // basedate // " " // basesec
    time = mdcur + mscur/secspday

    call clm_addvar(clmvar_double,ncid,'time',(/'time'/),'time',trim(str))

    ! Define surface grid (coordinate variables, latitude, longitude,
    ! surface type).

    call clm_addvar(clmvar_double,ncid,'longxy',(/'gridcell'/), &
          long_name='longitude', units='degrees_east')
    call clm_addvar(clmvar_double,ncid,'latixy',(/'gridcell'/), &
          long_name='latitude', units='degrees_north')
    call clm_addvar(clmvar_integer,ncid,'landmask',(/'gridcell'/), &
          long_name='land/ocean mask (0.=ocean and 1.=land)')

    ! Define time information

    call clm_addvar(clmvar_integer,ncid,'mcdate',(/'time'/), &
         long_name='current date (YYYYMMDD)')
    call clm_addvar(clmvar_integer,ncid,'mcsec',(/'time'/), &
         long_name='current seconds of current date', units='s')
    call clm_addvar(clmvar_integer,ncid,'mdcur',(/'time'/), &
         long_name='current day (from base day)')
    call clm_addvar(clmvar_integer,ncid,'mscur',(/'time'/), &
         long_name='current seconds of current day', units='s')
    call clm_addvar(clmvar_integer,ncid,'nstep',(/'time'/), &
         long_name='time step', units='s')

    ! Define time dependent variables

     call clm_addvar(clmvar_double,ncid,'FPCGRID', (/'pft ','time'/), &
          long_name='plant functional type cover', &
          units='fraction of vegetated area', &
          missing_value=1,fill_value=1)
     call clm_addvar(clmvar_double,ncid,'NIND', (/'pft ','time'/), &
          long_name='number of individuals', &
          units='individuals/m2 vegetated land', &
          missing_value=1,fill_value=1)

    call clm_enddef(ncid)

    ! -----------------------------------------------------------------------
    ! Write variables
    ! -----------------------------------------------------------------------

    ! Write surface grid (coordinate variables, latitude,
    ! longitude, surface type).

    call clm_writevar(ncid,'longxy',ldomain%lonc,gcomm_gridcell)
    call clm_writevar(ncid,'latixy',ldomain%latc,gcomm_gridcell)
    call clm_writevar(ncid,'landmask',ldomain%mask,gcomm_gridcell)

    ! Write current date, current seconds, current day, current nstep

    call curr_date(idatex,yr,mon,day,mcsec)
    mcdate = yr*10000 + mon*100 + day

    iktau = int(ktau,ik4)
    call clm_writevar(ncid,'mcdate',mcdate,nt=1)
    call clm_writevar(ncid,'mcsec',mcsec,nt=1)
    call clm_writevar(ncid,'mdcur',mdcur,nt=1)
    call clm_writevar(ncid,'mscur',mscur,nt=1)
    call clm_writevar(ncid,'nstep',iktau,nt=1)
    call clm_writevar(ncid,'time',time,nt=1)

    ! Write time dependent variables to CNDV history file

    ! The if .not. ifspecial statment below guarantees that the m index will
    ! always lie between 1 and maxpatch_pft

    rbuf2dg(:) = 0.D0
    do p = begp , endp
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(p) = fpcgrid(p)*100.D0
    end do
    call clm_writevar(ncid,'FPCGRID',rbuf2dg,gcomm_pft,nt=1)

    rbuf2dg(:) = 0.D0
    do p = begp , endp
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(p) = nind(p)
    end do
    call clm_writevar(ncid,'NIND',rbuf2dg,gcomm_pft,nt=1)

    ! Deallocate dynamic memory

    deallocate(rbuf2dg)

    !------------------------------------------------------------------
    ! Close and archive netcdf CNDV history file
    !------------------------------------------------------------------

    call clm_closefile(ncid)

    if (myid == italk) then
       write(stdout,*)'(histCNDV): Finished writing CNDV history dataset ',&
            trim(dgvm_fn), 'at nstep = ',ktau
    end if

  end subroutine histCNDV
  !
  ! Determine initial dataset filenames
  !
  character(len=256) function set_dgvm_filename ()
    implicit none
    character(len=256) :: cdate       !date char string
    integer(ik4) :: day                    !day (1 -> 31)
    integer(ik4) :: mon                    !month (1 -> 12)
    integer(ik4) :: yr                     !year (0 -> ...)
    integer(ik4) :: sec                    !seconds into current day
    call curr_date(idatex,yr,mon,day,sec)
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
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
    integer(ik4) , intent(in)  :: lbp , ubp
    ! number of pfts in naturally-vegetated filter
    integer(ik4) , intent(out) :: num_natvegp
    ! pft filter for naturally-vegetated points
    integer(ik4) , intent(out) :: filter_natvegp(ubp-lbp+1)
    logical , pointer :: present(:)     ! whether this pft present in patch
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
