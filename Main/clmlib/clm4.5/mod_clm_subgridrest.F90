module mod_clm_subgridrest
  use mod_intkinds
  use mod_realkinds
  use mod_clm_type
  use mod_mpmessage
  use mod_date
  use mod_runparams
  use mod_clm_nchelper
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
    integer(ik4) :: l , c , p ! indices
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
      call curr_date(idatex,yr,mon,day,mcsec)
      mcdate = yr*10000 + mon*100 + day
      call clm_writevar(ncid,'mcdate',mcdate)
      call clm_writevar(ncid,'mcsec',mcsec)
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
      call clm_writevar(ncid,'grid1d_lon',gptr%londeg)
      call clm_writevar(ncid,'grid1d_lat',gptr%latdeg)
    end if

    ! Write landunit info

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'land1d_lon', &
              (/'landunit'/),long_name='landunit longitude', &
              units='degrees_east')
      call clm_addvar(clmvar_double,ncid,'land1d_lat', &
              (/'landunit'/),long_name='landunit latitude', &
              units='degrees_north')
      call clm_addvar(clmvar_double,ncid,'land1d_wtxy',(/'landunit'/), &
              long_name='landunit weight relative to corresponding gridcell')
      call clm_addvar(clmvar_double,ncid,'land1d_ityplun',(/'landunit'/), &
           long_name='landunit type (vegetated,urban,lake,wetland or glacier)')
    else if (flag == 'write') then
      do l = begl , endl
        rlarr(l) = gptr%londeg(lptr%gridcell(l))
      end do
      call clm_writevar_par(ncid,'land1d_lon',rlarr,gcomm_landunit)
      do l = begl , endl
        rlarr(l) = gptr%latdeg(lptr%gridcell(l))
      end do
      call clm_writevar_par(ncid,'land1d_lat',rlarr,gcomm_landunit)
      call clm_writevar_par(ncid,'land1d_wtxy',lptr%wtgcell,gcomm_landunit)
      call clm_writevar_par(ncid,'land1d_ityplun',lptr%itype,gcomm_landunit)
    end if

    ! Write column info

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'cols1d_lon', &
              (/'column'/),long_name='column longitude', &
              units='degrees_east')
      call clm_addvar(clmvar_double,ncid,'cols1d_lat', &
              (/'column'/),long_name='column latitude', &
              units='degrees_north')
      call clm_addvar(clmvar_double,ncid,'cols1d_wtxy',(/'column'/), &
            long_name='column weight relative to corresponding gridcell')
      call clm_addvar(clmvar_double,ncid,'cols1d_wtlnd',(/'column'/), &
            long_name='column weight relative to corresponding landunit')
      call clm_addvar(clmvar_integer,ncid,'cols1d_ityplun',(/'column'/), &
           long_name='column landunit type (veget.,urban,lake,wetland,glacier)')
      call clm_addvar(clmvar_integer,ncid,'cols1d_ityp',(/'column'/), &
              long_name='column type (61-roof,62-sunwall,63-shadewall,'// &
                    '64-impervious road,65-pervious road,1-all other columns)')
    else if ( flag == 'write' ) then
      do c = begc , endc
        rcarr(c) = gptr%londeg(cptr%gridcell(c))
      end do
      call clm_writevar_par(ncid,'cols1d_lon',rcarr,gcomm_column)
      do c = begc , endc
        rcarr(c) = gptr%latdeg(cptr%gridcell(c))
      end do
      call clm_writevar_par(ncid,'cols1d_lat',rcarr,gcomm_column)
      call clm_writevar_par(ncid,'cols1d_wtxy',cptr%wtgcell,gcomm_column)
      call clm_writevar_par(ncid,'cols1d_wtlnd',cptr%wtlunit,gcomm_column)
      do c = begc , endc
        icarr(c) = lptr%itype(cptr%landunit(c))
      end do
      call clm_writevar_par(ncid,'cols1d_ityplun',icarr,gcomm_column)
      do c = begc , endc
        icarr(c) = cptr%itype((c))
      end do
      call clm_writevar_par(ncid,'cols1d_ityp',icarr,gcomm_column)
    end if

    ! Write pft info

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'pfts1d_lon', &
              (/'pft'/),long_name='pft longitude', &
              units='degrees_east')
      call clm_addvar(clmvar_double,ncid,'pfts1d_lat', &
              (/'pft'/),long_name='pft latitude', &
              units='degrees_north')
      call clm_addvar(clmvar_double,ncid,'pfts1d_wtxy', &
            (/'pft'/),long_name='pft weight relative to corresponding gridcell')
      call clm_addvar(clmvar_double,ncid,'pfts1d_wtlnd', &
            (/'pft'/),long_name='pft weight relative to corresponding landunit')
      call clm_addvar(clmvar_double,ncid,'pfts1d_wtcol', &
            (/'pft'/),long_name='pft weight relative to corresponding column')
      call clm_addvar(clmvar_integer,ncid,'pfts1d_itypveg', &
              (/'pft'/),long_name='pft vegetation type')
      call clm_addvar(clmvar_integer,ncid,'pfts1d_ityplun',(/'pft'/), &
          long_name='pft landunit type (vegetated,urban,lake,wetland,glacier)')
    else if ( flag == 'write' ) then
      do p = begp , endp
        rparr(p) = gptr%londeg(pptr%gridcell(p))
      end do
      call clm_writevar_par(ncid,'pfts1d_lon',rparr,gcomm_pft)
      do p = begp , endp
        rparr(p) = gptr%latdeg(pptr%gridcell(p))
      end do
      call clm_writevar_par(ncid,'pfts1d_lat',rparr,gcomm_pft)
      call clm_writevar_par(ncid,'pfts1d_wtxy',pptr%wtgcell,gcomm_pft)
      call clm_writevar_par(ncid,'pfts1d_wtlnd',pptr%wtlunit,gcomm_pft)
      call clm_writevar_par(ncid,'pfts1d_wtcol',pptr%wtcol,gcomm_pft)
      call clm_writevar_par(ncid,'pfts1d_itypveg',pptr%itype,gcomm_pft)
      do p = begp , endp
        iparr(p) = lptr%itype(pptr%landunit(p))
      end do
      call clm_writevar_par(ncid,'pfts1d_ityplun',iparr,gcomm_pft)
    end if

    if ( flag == 'write' ) then
      deallocate(rgarr,rlarr,rcarr,rparr)
      deallocate(igarr,ilarr,icarr,iparr)
    end if
  end subroutine subgridRest

end module mod_clm_subgridrest
