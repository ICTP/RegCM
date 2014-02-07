module subgridRestMod
  use mod_intkinds
  use mod_realkinds
  use mod_clm_type
  use mod_mpmessage
  use mod_date
  use mod_clm_domain , only : ldomain
  use mod_clm_decomp

  implicit none

  private

  public :: subgridRest

  contains

  subroutine subgridRest( ncid, flag )
    implicit none
    ! netCDF dataset id
    type(clm_filetype) , intent(inout) :: ncid
    ! flag to determine if define, write or read data
    character(len=*) , intent(in) :: flag
    integer(ik4) :: g , l , c , p , j , i ! indices
    integer(ik4) :: yr        ! current year (0 -> ...)
    integer(ik4) :: mon       ! current month (1 -> 12)
    integer(ik4) :: day       ! current day (1 -> 31)
    integer(ik4) :: mcsec     ! seconds of current date
    integer(ik4) :: mcdate    ! current date
    integer(ik4) :: begp , endp    ! per-proc beg/end pft indices
    integer(ik4) :: begc , endc    ! per-proc beg/end column indices
    integer(ik4) :: begl , endl    ! per-proc beg/end landunit indices
    integer(ik4) :: begg , endg    ! per-proc beg/end gridcell indices
    integer(ik4) :: ier            ! error status
    real(rk8) , pointer , dimension(:) :: rgarr     ! temporary
    real(rk8) , pointer , dimension(:) :: rlarr     ! temporary
    real(rk8) , pointer , dimension(:) :: rcarr     ! temporary
    real(rk8) , pointer , dimension(:) :: rparr     ! temporary
    integer(ik4) , pointer , dimension(:) :: igarr  ! temporary
    integer(ik4) , pointer , dimension(:) :: ilarr  ! temporary
    integer(ik4) , pointer , dimension(:) :: icarr  ! temporary
    integer(ik4) , pointer , dimension(:) :: iparr  ! temporary
    ! pointer to gridcell derived subtype
    type(gridcell_type) , pointer :: gptr
    ! pointer to landunit derived subtype
    type(landunit_type) , pointer :: lptr
    ! pointer to column derived subtype
    type(column_type) , pointer :: cptr
    ! pointer to pft derived subtype
    type(pft_type) , pointer :: pptr
    character(len=32) :: subname = 'SubgridRest' ! subroutine name

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Get relevant sizes

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! Allocate dynamic memory

    if ( flag == 'write') then
      allocate(rgarr(begg:endg), &
               rlarr(begl:endl), &
               rcarr(begc:endc), &
               rparr(begp:endp), stat=ier)
      if (ier /= 0) then
        call fatal(__FILE__,__LINE__, &
              'allocation error from inicfile_fields rarrs')
      end if
      allocate(igarr(begg:endg), &
               ilarr(begl:endl), &
               icarr(begc:endc), &
               iparr(begp:endp), stat=ier)
      if (ier /= 0) then
        call fatal(__FILE__,__LINE__, &
              'allocation error from inicfile_fields iarrs')
      end if
    end if

    ! Write output data (first write current date and seconds of current date)

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_integer,ncid,'mcdate', &
            long_name='current date as 8 digit integer (YYYYMMDD)')
      call clm_addvar(clmvar_integer,ncid,'mcsec', &
            long_name='current seconds of current date', units='s')
    else if ( flag == 'write' ) then
      call curr_date(yr,mon,day,mcsec)
      mcdate = yr*10000 + mon*100 + day
    end if

    ! Write gridcell info

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'grid1d_lon', &
              (/'gridcell'/),long_name='gridcell longitude', &
              units='degrees_east')
      call clm_addvar(clmvar_double,ncid,'grid1d_lat', &
              (/'gridcell'/),long_name='gridcell latitude', &
              units='degrees_north')
    else if ( flag == 'write' ) then
      call ncd_io(varname='grid1d_lon', data=gptr%londeg, dim1name=nameg, ncid=ncid, flag=flag)
      call ncd_io(varname='grid1d_lat', data=gptr%latdeg, dim1name=nameg, ncid=ncid, flag=flag)
    end if

    ! Write landunit info

    if (flag == 'define') then
      call ncd_defvar(ncid=ncid, varname='land1d_lon', xtype=ncd_double,  &
            dim1name='landunit', long_name='landunit longitude', units='degrees_east')
      call ncd_defvar(ncid=ncid, varname='land1d_lat', xtype=ncd_double,  &
            dim1name='landunit', long_name='landunit latitude', units='degrees_north')
      call ncd_defvar(ncid=ncid, varname='land1d_ixy', xtype=ncd_int,  &
            dim1name='landunit', long_name='2d longitude index of corresponding landunit')
      call ncd_defvar(ncid=ncid, varname='land1d_jxy', xtype=ncd_int,  &
            dim1name='landunit', long_name='2d latitude index of corresponding landunit')
      call ncd_defvar(ncid=ncid, varname='land1d_wtxy', xtype=ncd_double,  &
            dim1name='landunit', long_name='landunit weight relative to corresponding gridcell')
      call ncd_defvar(ncid=ncid, varname='land1d_ityplun', xtype=ncd_int,  &
            dim1name='landunit', long_name='landunit type (vegetated,urban,lake,wetland or glacier)')
    else if (flag == 'write') then
      do l=begl,endl
        rlarr(l) = gptr%londeg(lptr%gridcell(l))
      end do
      call ncd_io(varname='land1d_lon'    , data=rlarr        , dim1name=namel, ncid=ncid, flag=flag)
      do l=begl,endl
        rlarr(l) = gptr%latdeg(lptr%gridcell(l))
      end do
      call ncd_io(varname='land1d_lat'    , data=rlarr        , dim1name=namel, ncid=ncid, flag=flag)
      do l=begl,endl
        ilarr(l) = mod(ldecomp%gdc2glo(lptr%gridcell(l))-1,ldomain%ni) + 1
      end do
      call ncd_io(varname='land1d_ixy'    , data=ilarr        , dim1name=namel, ncid=ncid, flag=flag)
      do l=begl,endl
        ilarr(l) = (ldecomp%gdc2glo(lptr%gridcell(l))-1)/ldomain%ni + 1
      end do
      call ncd_io(varname='land1d_jxy'    , data=ilarr        , dim1name=namel, ncid=ncid, flag=flag)
      call ncd_io(varname='land1d_wtxy'   , data=lptr%wtgcell , dim1name=namel, ncid=ncid, flag=flag)
      call ncd_io(varname='land1d_ityplun', data=lptr%itype   , dim1name=namel, ncid=ncid, flag=flag)
    end if

    ! Write column info

    if (flag == 'define') then
      call ncd_defvar(ncid=ncid, varname='cols1d_lon', xtype=ncd_double,  &
           dim1name='column', long_name='column longitude', units='degrees_east')
      call ncd_defvar(ncid=ncid, varname='cols1d_lat', xtype=ncd_double,  &
           dim1name='column', long_name='column latitude', units='degrees_north')
      call ncd_defvar(ncid=ncid, varname='cols1d_ixy', xtype=ncd_int,   &
           dim1name='column', long_name='2d longitude index of corresponding column')
      call ncd_defvar(ncid=ncid, varname='cols1d_jxy', xtype=ncd_int,   &
           dim1name='column', long_name='2d latitude index of corresponding column')
      call ncd_defvar(ncid=ncid, varname='cols1d_wtxy', xtype=ncd_double,   &
            dim1name='column', long_name='column weight relative to corresponding gridcell')
      call ncd_defvar(ncid=ncid, varname='cols1d_wtlnd', xtype=ncd_double,   &
            dim1name='column', long_name='column weight relative to corresponding landunit')
      call ncd_defvar(ncid=ncid, varname='cols1d_ityplun', xtype=ncd_int,   &
            dim1name='column', long_name='column landunit type (vegetated,urban,lake,wetland or glacier)')
      call ncd_defvar(ncid=ncid, varname='cols1d_ityp', xtype=ncd_int,   &
            dim1name='column', long_name=&
           'column type (61-roof,62-sunwall,63-shadewall,64-impervious road,65-pervious road,1-all other columns)')
    else if (flag == 'write') then
      do c=begc,endc
        rcarr(c) = gptr%londeg(cptr%gridcell(c))
      end do
      call ncd_io(varname='cols1d_lon'  , data=rcarr        , dim1name=namec, ncid=ncid, flag=flag)
      do c=begc,endc
        rcarr(c) = gptr%latdeg(cptr%gridcell(c))
      end do
      call ncd_io(varname='cols1d_lat'  , data=rcarr        , dim1name=namec, ncid=ncid, flag=flag)
      do c=begc,endc
        icarr(c) = mod(ldecomp%gdc2glo(cptr%gridcell(c))-1,ldomain%ni) + 1
      end do
      call ncd_io(varname='cols1d_ixy'  , data=icarr        , dim1name=namec, ncid=ncid, flag=flag)
      do c=begc,endc
        icarr(c) = (ldecomp%gdc2glo(cptr%gridcell(c))-1)/ldomain%ni + 1
      end do
      call ncd_io(varname='cols1d_jxy'  , data=icarr        , dim1name=namec, ncid=ncid, flag=flag)
      call ncd_io(varname='cols1d_wtxy' , data=cptr%wtgcell , dim1name=namec, ncid=ncid, flag=flag)
      call ncd_io(varname='cols1d_wtlnd', data=cptr%wtlunit , dim1name=namec, ncid=ncid, flag=flag)
      do c=begc,endc
        icarr(c) = lptr%itype(cptr%landunit(c))
      end do
      call ncd_io(varname='cols1d_ityplun', data=icarr      , dim1name=namec, ncid=ncid, flag=flag)
      do c=begc,endc
        icarr(c) = cptr%itype((c))
      end do
      call ncd_io(varname='cols1d_ityp', data=icarr      , dim1name=namec, ncid=ncid, flag=flag)
    end if

    ! Write pft info

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pfts1d_lon', xtype=ncd_double,  &
            dim1name='pft', long_name='pft longitude', units='degrees_east')
       call ncd_defvar(ncid=ncid, varname='pfts1d_lat', xtype=ncd_double,  &
            dim1name='pft', long_name='pft latitude', units='degrees_north')
       call ncd_defvar(ncid=ncid, varname='pfts1d_ixy', xtype=ncd_int,  &
            dim1name='pft', long_name='2d longitude index of corresponding pft')
       call ncd_defvar(ncid=ncid, varname='pfts1d_jxy', xtype=ncd_int,  &
            dim1name='pft', long_name='2d latitude index of corresponding pft')
       call ncd_defvar(ncid=ncid, varname='pfts1d_ci', xtype=ncd_int,  &
            dim1name='pft', long_name='1d column index of corresponding pft')
       call ncd_defvar(ncid=ncid, varname='pfts1d_wtxy', xtype=ncd_double,  &
            dim1name='pft', long_name='pft weight relative to corresponding gridcell')
       call ncd_defvar(ncid=ncid, varname='pfts1d_wtlnd', xtype=ncd_double,  &
            dim1name='pft', long_name='pft weight relative to corresponding landunit')
       call ncd_defvar(ncid=ncid, varname='pfts1d_wtcol', xtype=ncd_double,  &
            dim1name='pft', long_name='pft weight relative to corresponding column')
       call ncd_defvar(ncid=ncid, varname='pfts1d_itypveg', xtype=ncd_int,  &
            dim1name='pft', long_name='pft vegetation type')
       call ncd_defvar(ncid=ncid, varname='pfts1d_ityplun', xtype=ncd_int,  &
            dim1name='pft', long_name='pft landunit type (vegetated,urban,lake,wetland or glacier)')
    else if (flag == 'write') then
       do p=begp,endp
          rparr(p) = gptr%londeg(pptr%gridcell(p))
       end do
       call ncd_io(varname='pfts1d_lon'    , data=rparr        , dim1name=namep, ncid=ncid, flag=flag)
       do p=begp,endp
          rparr(p) = gptr%latdeg(pptr%gridcell(p))
       end do
       call ncd_io(varname='pfts1d_lat'    , data=rparr        , dim1name=namep, ncid=ncid, flag=flag)
      do p=begp,endp
        iparr(p) = mod(ldecomp%gdc2glo(pptr%gridcell(p))-1,ldomain%ni) + 1
      end do
      call ncd_io(varname='pfts1d_ixy'    , data=iparr        , dim1name=namep, ncid=ncid, flag=flag)
      do p=begp,endp
        iparr(p) = (ldecomp%gdc2glo(pptr%gridcell(p))-1)/ldomain%ni + 1
      end do
      call ncd_io(varname='pfts1d_jxy'    , data=iparr        , dim1name=namep, ncid=ncid, flag=flag)
      call ncd_io(varname='pfts1d_wtxy'   , data=pptr%wtgcell , dim1name=namep, ncid=ncid, flag=flag)
      call ncd_io(varname='pfts1d_wtlnd'  , data=pptr%wtlunit , dim1name=namep, ncid=ncid, flag=flag)
      call ncd_io(varname='pfts1d_wtcol'  , data=pptr%wtcol   , dim1name=namep, ncid=ncid, flag=flag)
      call ncd_io(varname='pfts1d_itypveg', data=pptr%itype   , dim1name=namep, ncid=ncid, flag=flag)
      do p=begp,endp
        iparr(p) = lptr%itype(pptr%landunit(p))
      end do
      call ncd_io(varname='pfts1d_ityplun', data=iparr      , dim1name=namep, ncid=ncid, flag=flag)
    end if

    if ( flag == 'write' ) then
      deallocate(rgarr,rlarr,rcarr,rparr)
      deallocate(igarr,ilarr,icarr,iparr)
    end if
  end subroutine subgridRest

end module subgridRestMod
