module mod_clm_cndv

#if (defined CNDV)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNDVMod
!
! !DESCRIPTION:
! Module containing routines to drive the annual dynamic vegetation
! that works with CN, reset related variables,
! and initialize/reset time invariant variables
!
! !USES:
  use mod_realkinds
  use mod_dynparam
  use mod_mppparam
  use mod_clm_cnvegstructupdate, only : CNVegStructUpdate
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public dv                 ! Drives the annual dynamic vegetation that
                            ! works with CN
  public histCNDV           ! Output CNDV history file
!
! !REVISION HISTORY:
! Module modified by Sam Levis from similar module DGVMMod
! created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: dv
!
! !INTERFACE:
  subroutine dv(lbg, ubg, lbp, ubp, num_natvegp, filter_natvegp, kyr)
!
! !DESCRIPTION:
! Drives the annual dynamic vegetation that works with CN
!
! !USES:
    use mod_clm_type
    use mod_clm_cndvlight        , only : Light
    use mod_clm_cndvestablishment, only : Establishment
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbg, ubg       ! gridcell bounds
    integer, intent(in) :: lbp, ubp       ! pft bounds
    integer, intent(inout) :: num_natvegp ! number of naturally-vegetated
                                          ! pfts in filter
    integer, intent(inout) :: filter_natvegp(ubp-lbp+1) ! filter for
                                          ! naturally-vegetated pfts
    integer, intent(in) :: kyr            ! used in routine climate20 below
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
   integer , pointer :: mxy(:)         ! pft m index (for laixy(i,j,m),etc.)
   integer , pointer :: pgridcell(:)   ! gridcell of corresponding pft
   real(rk8), pointer :: fpcgrid(:)     ! foliar projective cover on gridcell (fraction)
   real(rk8), pointer :: agdd(:)        ! accumulated growing degree days above 5
   real(rk8), pointer :: t_mo_min(:)    ! annual min of t_mo (Kelvin)
!
! local pointers to implicit inout arguments
!
   real(rk8), pointer :: tmomin20(:)         ! 20-yr running mean of tmomin
   real(rk8), pointer :: agdd20(:)           ! 20-yr running mean of agdd
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: g,p                    ! indices
!-----------------------------------------------------------------------

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

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: histCNDV
!
! !INTERFACE:
  subroutine histCNDV()
!
! !DESCRIPTION:
! Create CNDV history file
!
! !USES:
    use mod_clm_clmtype
    use mod_clm_decomp       , only : get_proc_bounds, get_proc_global
    use mod_clm_varpar      , only : maxpatch_pft
    use mod_clm_domain       , only : ldomain
    use mod_clm_varctl      , only : caseid, ctitle, finidat, fsurdat, fpftcon, iulog
    use mod_clm_varcon      , only : spval
    use mod_clm_time_manager, only : get_ref_date, get_nstep, get_curr_date, get_curr_time
    use mod_clm_varcon  , only : secspday
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
   logical , pointer :: ifspecial(:)        ! true=>landunit is not vegetated (landunit-level)
   integer , pointer :: pgridcell(:)        ! gridcell index of corresponding pft (pft-level)
   integer , pointer :: plandunit(:)        ! landunit index of corresponding pft (pft-level)
   integer , pointer :: mxy(:)              ! pft m index (for laixy(i,j,m),etc.)
   real(rk8), pointer :: fpcgrid(:)          ! foliar projective cover on gridcell (fraction)
   real(rk8), pointer :: nind(:)             ! number of individuals (#/m**2)
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: dgvm_fn      ! dgvm history filename
    integer :: ncid                    ! netcdf file id
    integer :: g,p,l                   ! indices
    integer :: begp, endp              ! per-proc beginning and ending pft indices
    integer :: begc, endc              ! per-proc beginning and ending column indices
    integer :: begl, endl              ! per-proc beginning and ending landunit indices
    integer :: begg, endg              ! per-proc gridcell ending gridcell indices
    integer :: ihost
    integer :: ier                     ! error status
    integer :: mdcur, mscur, mcdate    ! outputs from get_curr_time
    integer :: yr,mon,day,mcsec        ! outputs from get_curr_date
    integer :: hours,minutes,secs      ! hours,minutes,seconds of hh:mm:ss
    integer :: nstep                   ! time step
    integer :: nbsec                   ! seconds components of a date
    integer , dimension(5) :: dimids   ! dimension, variable id
    integer , dimension(5) :: varids   ! variable id
    integer :: tid
    integer :: ipnt
    real(rk8):: time                   ! current time
    character(len=256) :: str          ! temporary string
    character(len=  8) :: curdate      ! current date
    character(len=  8) :: curtime      ! current time
    character(len= 10) :: basedate     ! base date (yyyymmdd)
    character(len=  8) :: basesec      ! base seconds
    real(rk8), pointer :: rbuf2dg(:,:)  ! temporary
    character(len=32) :: subname='histCNDV'
    character (len=32) :: hostname='?'
    character (len=32) :: user='?'
!-----------------------------------------------------------------------

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

    allocate(rbuf2dg(begg:endg,maxpatch_pft), stat=ier)
    if (ier /= 0) call fatal(__FILE__,__LINE__,&
      'histCNDV: allocation error for rbuf2dg')

    ! -----------------------------------------------------------------------
    ! Create new netCDF file. File will be in define mode
    ! -----------------------------------------------------------------------

    dgvm_fn = set_dgvm_filename()

    call createfile_withname(trim(dgvm_fn),ncid)

    ! -----------------------------------------------------------------------
    ! Create global attributes.
    ! -----------------------------------------------------------------------
    
    str = 'CF1.0'
    call add_attribute(ncid, 'conventions', trim(str))
    
    call getdatetime(curdate, curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call add_attribute(ncid, 'history', trim(str))
 
#ifdef IBM
    hostname='ibm platform '
    user= 'Unknown'
#else
    ihost = hostnm(hostname)
    call getlog(user)
#endif
       
    call add_attribute(ncid, 'logname', trim(user))
    call add_attribute(ncid, 'host', trim(hostname))
       
    str = 'Community Land Model: CLM3'
    call add_attribute(ncid, 'source',  trim(str))
       
    str = '$Name$'
    call add_attribute(ncid, 'version', trim(str))
       
    str = '$Id$'
    call add_attribute(ncid, 'revision_id',  trim(str))

    str = ctitle
    call add_attribute(ncid, 'case_title', trim(str))

    str = caseid
    call add_attribute(ncid, 'case_id',  trim(str))

    call add_attribute(ncid, 'Surface_dataset',  trim(fsurdat))

    str = 'arbitrary initialization'
    if (finidat /= ' ') str = finidat
    call add_attribute(ncid, 'Initial_conditions_dataset',  trim(str))

    call add_attribute(ncid,'PFT_physiological_constants_dataset',trim(fpftcon))

    ! -----------------------------------------------------------------------
    ! Define dimensions.
    ! -----------------------------------------------------------------------
    
    ipnt = 1
    if (ldomain%isgrid2d) then
      call add_dimension(ncid,'lon',ldomain%ni,ipnt,dimids)
      call add_dimension(ncid,'lat',ldomain%nj,ipnt,dimids)
    else
      call add_dimension(ncid,'gridcell',ldomain%ns,ipnt,dimids)
    end if
    call add_dimension(ncid,'pft',maxpatch_pft,ipnt,dimids)
    tid = ipnt
    call add_dimension(ncid,'time',-1,ipnt,dimids)
    call add_dimension(ncid,'string_length',80,ipnt,dimids)
    
    ! -----------------------------------------------------------------------
    ! Define variables
    ! -----------------------------------------------------------------------
    
    ! Define coordinate variables (including time)
    ipnt = 1
    if (ldomain%isgrid2d) then
       call add_variable(ncid,'lon','coordinate longitude','degrees_east', &
         dimids(1:1),ipnt,varids)
       call add_variable(ncid,'lat','coordinate latitude','degrees_north', &
         dimids(2:2),ipnt,varids)
    end if
    
    call get_curr_time(mdcur, mscur)
    call get_ref_date(yr, mon, day, nbsec)
    hours   = nbsec / 3600
    minutes = (nbsec - hours*3600) / 60
    secs    = (nbsec - hours*3600 - minutes*60)
    write(basedate,80) yr,mon,day
80  format(i4.4,'-',i2.2,'-',i2.2)
    write(basesec ,90) hours, minutes, secs
90  format(i2.2,':',i2.2,':',i2.2)
    str = 'days since ' // basedate // " " // basesec
    time = mdcur + mscur/secspday
    
    call add_variable(ncid,'time','time',trim(str),dimids(tid:tid),ipnt,varids)
       
    ! Define surface grid (coordinate variables, latitude, longitude, surface type).
    
    if (ldomain%isgrid2d) then
       call ncd_defvar(ncid=ncid, varname='longxy', xtype=ncprec, &
            dim1name='lon', dim2name='lat', &
            long_name='longitude', units='degrees_east')
       
       call ncd_defvar(ncid=ncid, varname='latixy', xtype=ncprec, &
            dim1name='lon', dim2name='lat', &
            long_name='latitude', units='degrees_north')
       
       call ncd_defvar(ncid=ncid, varname='landmask', xtype=ncd_int, &
            dim1name='lon', dim2name='lat', &
            long_name='land/ocean mask (0.=ocean and 1.=land)')
    else
       call ncd_defvar(ncid=ncid, varname='longxy', xtype=ncprec, &
            dim1name='gridcell',&
            long_name='longitude', units='degrees_east')
       
       call ncd_defvar(ncid=ncid, varname='latixy', xtype=ncprec, &
            dim1name='gridcell',&
            long_name='latitude', units='degrees_north')
       
       call ncd_defvar(ncid=ncid, varname='landmask', xtype=ncd_int, &
            dim1name='gridcell', &
            long_name='land/ocean mask (0.=ocean and 1.=land)')
    end if

    ! Define time information

    call ncd_defvar(ncid=ncid, varname='mcdate', xtype=ncd_int, dim1name='time',&
         long_name='current date (YYYYMMDD)')

    call ncd_defvar(ncid=ncid, varname='mcsec', xtype=ncd_int, dim1name='time',&
         long_name='current seconds of current date', units='s')

    call ncd_defvar(ncid=ncid, varname='mdcur', xtype=ncd_int, dim1name='time',&
         long_name='current day (from base day)')

    call ncd_defvar(ncid=ncid, varname='mscur', xtype=ncd_int, dim1name='time',&
         long_name='current seconds of current day', units='s')

    call ncd_defvar(ncid=ncid, varname='nstep', xtype=ncd_int, dim1name='time',&
         long_name='time step', units='s')

    ! Define time dependent variables

    if (ldomain%isgrid2d) then
       call ncd_defvar(ncid=ncid, varname='FPCGRID', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='plant functional type cover', units='fraction of vegetated area', &
            missing_value=spval, fill_value=spval)
       
       call ncd_defvar(ncid=ncid, varname='NIND', xtype=ncprec, &
            dim1name='lon', dim2name='lat', dim3name='pft', dim4name='time', &
            long_name='number of individuals', units='individuals/m2 vegetated land', &
            missing_value=spval, fill_value=spval)
    else 
       call ncd_defvar(ncid=ncid, varname='FPCGRID', xtype=ncprec, &
            dim1name='gridcell', dim2name='pft', dim3name='time', &
            long_name='plant functional type cover', units='fraction of vegetated area', &
            missing_value=spval, fill_value=spval)
       
       call ncd_defvar(ncid=ncid, varname='NIND', xtype=ncprec, &
            dim1name='gridcell', dim2name='pft', dim3name='time', &
            long_name='number of individuals', units='individuals/m2 vegetated land', &
            missing_value=spval, fill_value=spval)
    end if

    call ncd_enddef(ncid)

    ! -----------------------------------------------------------------------
    ! Write variables
    ! -----------------------------------------------------------------------

    ! Write surface grid (coordinate variables, latitude, longitude, surface type).

    call ncd_io(ncid=ncid, varname='longxy'  , data=ldomain%lonc, flag='write', &
         dim1name=grlnd)
    call ncd_io(ncid=ncid, varname='latixy'  , data=ldomain%latc, flag='write', &
         dim1name=grlnd)
    call ncd_io(ncid=ncid, varname='landmask', data=ldomain%mask, flag='write', &
         dim1name=grlnd)

    ! Write current date, current seconds, current day, current nstep

    call get_curr_date(yr, mon, day, mcsec)
    mcdate = yr*10000 + mon*100 + day
    nstep = get_nstep()

    call ncd_io(ncid=ncid, varname='mcdate', data=mcdate, nt=1, flag='write')
    call ncd_io(ncid=ncid, varname='mcsec' , data=mcsec , nt=1, flag='write')
    call ncd_io(ncid=ncid, varname='mdcur' , data=mdcur , nt=1, flag='write')
    call ncd_io(ncid=ncid, varname='mscur' , data=mcsec , nt=1, flag='write')
    call ncd_io(ncid=ncid, varname='nstep' , data=nstep , nt=1, flag='write')
    call ncd_io(ncid=ncid, varname='time'  , data=time  , nt=1, flag='write')

    ! Write time dependent variables to CNDV history file

    ! The if .not. ifspecial statment below guarantees that the m index will
    ! always lie between 1 and maxpatch_pft

    rbuf2dg(:,:) = 0.D0
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(g,mxy(p)) = fpcgrid(p)*100.D0
    end do
    call ncd_io(ncid=ncid, varname='FPCGRID', dim1name=grlnd, data=rbuf2dg, &
         nt=1, flag='write')

    rbuf2dg(:,:) = 0.D0
    do p = begp,endp
       g = pgridcell(p)
       l = plandunit(p)
       if (.not. ifspecial(l)) rbuf2dg(g,mxy(p)) = nind(p)
    end do
    call ncd_io(ncid=ncid, varname='NIND', dim1name=grlnd, data=rbuf2dg, &
         nt=1, flag='write')

    ! Deallocate dynamic memory

    deallocate(rbuf2dg)

    !------------------------------------------------------------------
    ! Close and archive netcdf CNDV history file
    !------------------------------------------------------------------

    call closefile(ncid)

    if (myid == italk) then
       write(iulog,*)'(histCNDV): Finished writing CNDV history dataset ',&
            trim(dgvm_fn), 'at nstep = ',get_nstep()
    end if

  end subroutine histCNDV

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_dgvm_filename
!
! !INTERFACE:
  character(len=256) function set_dgvm_filename ()
!
! !DESCRIPTION:
! Determine initial dataset filenames
!
! !USES:
    use clm_varctl       , only : caseid, inst_suffix
    use clm_time_manager , only : get_curr_date
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
!-----------------------------------------------------------------------

    call get_curr_date (yr, mon, day, sec)
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    set_dgvm_filename = "./"//trim(caseid)//".clm2"//trim(inst_suffix)//&
                        ".hv."//trim(cdate)//".nc"

  end function set_dgvm_filename

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BuildNatVegFilter
!
! !INTERFACE:
  subroutine BuildNatVegFilter(lbp, ubp, num_natvegp, filter_natvegp)
!
! !DESCRIPTION:
! Reconstruct a filter of naturally-vegetated PFTs for use in DGVM
!
! !USES:
    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: lbp, ubp                   ! pft bounds
    integer, intent(out) :: num_natvegp                ! number of pfts in naturally-vegetated filter
    integer, intent(out) :: filter_natvegp(ubp-lbp+1)  ! pft filter for naturally-vegetated points
!
! !CALLED FROM:
! subroutine lpj in this module
!
! !REVISION HISTORY:
! Author: Forrest Hoffman
!
! !LOCAL VARIABLES:
! local pointers to implicit in arguments
    logical , pointer :: present(:)     ! whether this pft present in patch
!EOP
!
! !LOCAL VARIABLES:
    integer :: p
!-----------------------------------------------------------------------

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
