module mod_clm_biogeophysrest
  !
  ! Reads from or biogeophysics restart/initial data
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_stdio
  use mod_dynparam
  use mod_runparams
  use mod_mppparam
  use mod_clm_nchelper
  use mod_clm_type
  use mod_clm_decomp , only : get_proc_bounds , gcomm_gridcell , &
          gcomm_landunit , gcomm_column , gcomm_pft
  use mod_clm_varpar , only : nlevgrnd, nlevsno, nlevlak, nlevurb, nlevsoi, &
          nlevcan
  use mod_clm_varcon , only : istcrop
  use mod_clm_varcon , only : denice , denh2o , istdlak , istslak , isturb , &
          istsoil , pondmx , watmin , spval , icol_roof , icol_sunwall,      &
          icol_shadewall
  use mod_clm_varctl , only : allocate_all_vegpfts, nsrest, &
          pertlim , nsrContinue , nsrStartup
  use mod_clm_initsurfalb , only : do_initsurfalb
  use mod_clm_snicar , only : snw_rds_min
  use mod_clm_mkarbinit , only : perturbIC
  use mod_clm_atmlnd , only : clm_a2l

  implicit none

  private

  save

  public :: BiogeophysRest

  private :: weights_exactly_the_same
  private :: weights_within_roundoff_different
  private :: weights_tooDifferent

  contains
  !
  ! Read/Write biogeophysics information to/from restart file.
  !
  subroutine BiogeophysRest( ncid, flag )
    implicit none
    type(clm_filetype) , intent(inout) :: ncid ! netcdf id
    character(len=*) , intent(in)    :: flag ! 'read' or 'write'
    real(rk8) :: maxwatsat                 !maximum porosity
    real(rk8) :: excess                    !excess volumetric soil water
    real(rk8) :: totwat                    !total soil water (mm)
    real(rk8) :: maxdiff                   !maximum difference in PFT weights
    real(rk8) , pointer :: wtgcell(:)      ! Grid cell weights for PFT
    real(rk8) , pointer :: wtlunit(:)      ! Land-unit weights for PFT
    real(rk8) , pointer :: wtcol(:)        ! Column weights for PFT
    integer(ik4) :: p , c , l , j , iv ! indices
    real(rk8) , pointer :: zi(:,:)    ! interface level below a "z" level (m)
    integer(ik4) :: nlevs       ! number of layers
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    logical :: readvar      ! determine if variable is on initial file
    integer(ik4) , pointer :: clandunit(:) ! landunit of corresponding column
    integer(ik4) , pointer :: ltype(:)     ! landunit type
    integer(ik4) , pointer :: ctype(:)     ! column type
    type(landunit_type), pointer :: lptr ! pointer to landunit derived subtype
    type(column_type) , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type) , pointer :: pptr     ! pointer to pft derived subtype
    real(rk8) , pointer   :: temp2d(:,:) ! temporary for zisno
    ! tolerance of acceptible difference
    real(rk8) , parameter :: adiff = 5.e-4_rk8
    character(len=7) :: filetypes(0:3)
    character(len=32) :: fileusing
    character(len=*) , parameter :: sub = "BiogeophysRest"
    logical :: lstart

    filetypes(:)           = "missing"
    filetypes(nsrStartup)  = "finidat"
    filetypes(nsrContinue) = "restart"

    lstart = rcmtimer%integrating( )

    ! Set pointers into derived type

    lptr       => clm3%g%l
    cptr       => clm3%g%l%c
    pptr       => clm3%g%l%c%p
    ltype      => lptr%itype
    clandunit  => cptr%landunit
    clandunit  => cptr%landunit
    zi         => clm3%g%l%c%cps%zi
    ctype      => cptr%itype

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    !
    ! Read in weights if allocating all vegetation types
    !

    if ( allocate_all_vegpfts ) then

      ! pft weight wrt gridcell

      if ( flag == 'define' ) then
        call clm_addvar(clmvar_double,ncid,'PFT_WTGCELL',['pft'], &
              long_name='pft weight relative to corresponding gridcell', &
              units='')
      else if ( flag == 'read' ) then
        if ( lstart .and. .not. clm_check_var(ncid,'PFT_WTGCELL') ) then
          call fatal(__FILE__,__LINE__,'clm_now_stopping')
        else
          ! Copy weights calculated from fsurdat/fpftdyn to temp array
          ! for comparision
          ! Don't read directly into temp array -- so that answers are
          ! identical with clm3.6.58. EBK 1/9/2010
          allocate( wtgcell(begp:endp) )
          wtgcell(:) = pptr%wtgcell(:)
          call clm_readvar(ncid,'PFT_WTGCELL',pptr%wtgcell,gcomm_pft)
        end if
      else if (flag == 'write' ) then
        call clm_writevar(ncid,'PFT_WTGCELL',pptr%wtgcell,gcomm_pft)
      end if

      ! pft weight wrt landunit

      if ( flag == 'define' ) then
        call clm_addvar(clmvar_double,ncid,'PFT_WTLUNIT',['pft'], &
              long_name='pft weight relative to corresponding landunit', &
              units='')
      else if ( flag == 'read' ) then
        if ( lstart .and. .not. clm_check_var(ncid,'PFT_WTLUNIT') ) then
          call fatal(__FILE__,__LINE__,'clm_now_stopping')
        else
          ! Copy weights calculated from fsurdat/fpftdyn to temp array
          ! for comparision
          ! Don't read directly into temp array -- so that answers are
          ! identical with clm3.6.58. EBK 1/9/2010
          allocate( wtlunit(begp:endp) )
          wtlunit(:) = pptr%wtlunit(:)
          call clm_readvar(ncid,'PFT_WTLUNIT',pptr%wtlunit,gcomm_pft)
        end if
      else if ( flag == 'write' ) then
        call clm_writevar(ncid,'PFT_WTLUNIT',pptr%wtlunit,gcomm_pft)
      end if

      ! pft weight wrt column

      if ( flag == 'define' ) then
        call clm_addvar(clmvar_double,ncid,'PFT_WTCOL',['pft'], &
               long_name='pft weight relative to corresponding column', &
               units='')
      else if ( flag == 'read' ) then
        if ( lstart .and. .not. clm_check_var(ncid,'PFT_WTCOL') ) then
          call fatal(__FILE__,__LINE__,'clm_now_stopping')
        else
          ! Copy weights calculated from fsurdat/fpftdyn to temp array
          ! for comparision
          ! Don't read directly into temp array -- so that answers are
          ! identical with clm3.6.58. EBK 1/9/2010
          allocate( wtcol(begp:endp)   )
          wtcol(:) = pptr%wtcol(:)
          call clm_readvar(ncid,'PFT_WTCOL',pptr%wtcol,gcomm_pft)
        end if
      else if ( flag == 'write' ) then
        call clm_writevar(ncid,'PFT_WTCOL',pptr%wtcol,gcomm_pft)
      end if

      if ( flag == 'read' ) then
#ifdef DYNPFT
        fileusing = "fsurdat/fpftdyn"
#else
        fileusing = "fsurdat"
#endif
        !
        ! Note: Do not compare weights if restart
        !
        if ( nsrest == nsrContinue .or. fileusing == "fsurdat/fpftdyn" ) then
          ! Do NOT do any testing for restart or a pftdyn case
          !
          ! Otherwise test and make sure weights agree to reasonable tolerence
          !
        else if ( .not. weights_exactly_the_same(pptr,wtgcell, &
                                                 wtlunit,wtcol) ) then
#if (!defined CNDV)

          if ( weights_within_roundoff_different(pptr,wtgcell, &
                                                 wtlunit,wtcol) ) then
            write(stderr,*) &
               sub//"::NOTE, PFT weights from ", filetypes(nsrest),      &
                    " file and ", trim(fileusing), &
                    " file(s) are different to roundoff -- using ", &
                    trim(fileusing), " values."
          else if (weights_tooDifferent(begp,endp,pptr, &
                                        wtgcell,adiff,maxdiff) ) then
            write(stderr,*) &
            "ERROR:: PFT weights are SIGNIFICANTLY different from the input ", &
            filetypes(nsrest), " file and ", trim(fileusing), " file(s)."
            write(stderr,*) &
            "ERROR:: maximum difference is ", maxdiff, " max allowed = ", adiff
            write(stderr,*) &
            "ERROR:: Run interpinic on your initial condition file to "// &
            "interpolate to the new surface dataset"
            call fatal(__FILE__,__LINE__,&
                  sub//"::ERROR:: Weights between initial condition "// &
                       "file and surface dataset are too different")
          else
            write(stderr,*) sub//"::NOTE, PFT weights from ", &
                    filetypes(nsrest), " file and ", trim(fileusing), &
                   " file(s) are different to < ", &
                   adiff, " -- using ", trim(fileusing), " values."
          end if
          write(stderr,*) &
            sub//"::WARNING, weights different between ", filetypes(nsrest), &
            " file and ", trim(fileusing), &
            " file(s), but close enough -- using ",    &
            trim(fileusing), " values."
          ! Copy weights from fsurdat file back in -- they are only off
          ! by roundoff to 1% or so...
          pptr%wtgcell(:) = wtgcell(:)
          pptr%wtlunit(:) = wtlunit(:)
          pptr%wtcol(:)   = wtcol(:)
#endif
        end if
        deallocate( wtgcell )
        deallocate( wtlunit )
        deallocate( wtcol   )
      end if
    end if

    ! Note - for the snow interfaces, are only examing the snow interfaces
    ! above zi=0 which is why zisno and zsno have the same level dimension below
    ! (Note - for zisno, zi(0) is set to 0 in routine iniTimeConst)

    ! pft energy flux - eflx_lwrad_out

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'EFLX_LWRAD_OUT',['pft'], &
           long_name='emitted infrared (longwave) radiation', units='watt/m^2')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'EFLX_LWRAD_OUT') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'EFLX_LWRAD_OUT', &
                pptr%pef%eflx_lwrad_out,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'EFLX_LWRAD_OUT', &
               pptr%pef%eflx_lwrad_out,gcomm_pft)
    end if

    ! column water state variable - snow levels

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_integer,ncid,'SNLSNO',['column'], &
            long_name='number of snow layers', units='unitless')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'SNLSNO') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'SNLSNO',cptr%cps%snl,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'SNLSNO',cptr%cps%snl,gcomm_column)
    end if

    ! column water state variable - snow_depth
    ! As of clm4_0_76, SNOW_DEPTH is written to restarts.

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'SNOW_DEPTH',['column'], &
            long_name='snow depth', units='m')
    else if ( flag == 'read' ) then
      if ( lstart  .and. .not. clm_check_var(ncid,'SNOW_DEPTH') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'SNOW_DEPTH',cptr%cps%snow_depth,gcomm_column)
      end if
    else if ( flag == 'write' ) then
       call clm_writevar(ncid,'SNOW_DEPTH',cptr%cps%snow_depth,gcomm_column)
    end if

    ! column water state variable - int_snow

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'INT_SNOW',['column'], &
            long_name='accumulated snow', units='mm')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'INT_SNOW') ) then
        if ( rcmtimer%start( ) ) then
          cptr%cws%int_snow(:) = 0.0_rk8
        else
          call fatal(__FILE__,__LINE__,'clm_now_stopping')
        end if
      else
        call clm_readvar(ncid,'INT_SNOW',cptr%cws%int_snow,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'INT_SNOW',cptr%cws%int_snow,gcomm_column)
    end if

    ! column water state variable - wa

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'WA',['column'], &
            long_name='water in the unconfined aquifer', units='mm')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'WA') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'WA',cptr%cws%wa,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'WA',cptr%cws%wa,gcomm_column)
    end if

    ! gridcell type water flux variable - tws

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'TWS',['gridcell'], &
             long_name='total water storage', units='mm/s')
    else if ( flag == 'read' ) then
      readvar = clm_check_var(ncid,'TWS')
      if ( lstart .and. .not. readvar ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        if ( readvar ) then
          call clm_readvar(ncid,'TWS',clm3%g%tws,gcomm_gridcell)
        else
          ! initial run, not restart: initialize flood to zero
          clm3%g%tws = 0._rk8
        end if
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'TWS',clm3%g%tws,gcomm_gridcell)
    end if

    ! column water state variable - zwt

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'ZWT',['column'], &
            long_name='water table depth', units='m')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'ZWT') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'ZWT',cptr%cws%zwt,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'ZWT',cptr%cws%zwt,gcomm_column)
    end if

    ! column water state variable - frost_table

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'FROST_TABLE',['column'], &
            long_name='frost table depth', units='m')
    else if ( flag == 'read' ) then
      readvar = clm_check_var(ncid,'FROST_TABLE')
      if ( lstart .and. .not. readvar ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        if ( readvar ) then
          call clm_readvar(ncid,'FROST_TABLE',cptr%cws%frost_table,gcomm_column)
        else
          do c = begc , endc
            cptr%cws%frost_table(c) = zi(c,nlevsoi)
          end do
        end if
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'FROST_TABLE',cptr%cws%frost_table,gcomm_column)
    end if

    ! column water state variable - zwt_perched

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'ZWT_PERCH',['column'], &
            long_name='perched water table depth', units='m')
    else if ( flag == 'read' ) then
      readvar = clm_check_var(ncid,'ZWT_PERCH')
      if ( lstart .and. .not. readvar ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        if ( readvar ) then
          call clm_readvar(ncid,'ZWT_PERCH',cptr%cws%zwt_perched,gcomm_column)
        else
          do c = begc , endc
            cptr%cws%zwt_perched(c) = zi(c,nlevsoi)
          end do
        end if
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'ZWT_PERCH',cptr%cws%zwt_perched,gcomm_column)
    end if

    ! column type physical state variable - frac_sno_eff

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'frac_sno_eff',['column'], &
            long_name='fraction of ground covered by snow (0 to 1)', &
            units='unitless')
    else if ( flag == 'read' ) then
      readvar = clm_check_var(ncid,'frac_sno_eff')
      if ( lstart .and. .not. readvar ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        if ( readvar ) then
          call clm_readvar(ncid,'frac_sno_eff', &
                  cptr%cps%frac_sno_eff,gcomm_column)
        else
          do c = begc , endc
            cptr%cps%frac_sno_eff(c) = 0.0_rk8
          end do
        end if
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'frac_sno_eff',cptr%cps%frac_sno_eff,gcomm_column)
    end if

    ! column type physical state variable - frac_sno

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'frac_sno',['column'], &
            long_name='fraction of ground covered by snow (0 to 1)', &
            units='unitless')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'frac_sno') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'frac_sno',cptr%cps%frac_sno,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'frac_sno',cptr%cps%frac_sno,gcomm_column)
    end if

    ! column type physical state variable - dzsno

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'DZSNO',['column','levsno'], &
            long_name='snow layer thickness', units='m', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'DZSNO') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        allocate(temp2d(begc:endc,-nlevsno+1:0))
        call clm_readvar(ncid,'DZSNO',temp2d,gcomm_column, switchdim=.true.)
        cptr%cps%dz(begc:endc,-nlevsno+1:0) = temp2d(begc:endc,-nlevsno+1:0)
        deallocate(temp2d)
      end if
    else if ( flag == 'write' ) then
      allocate(temp2d(begc:endc,-nlevsno+1:0))
      temp2d(begc:endc,-nlevsno+1:0) = cptr%cps%dz(begc:endc,-nlevsno+1:0)
      call clm_writevar(ncid,'DZSNO',temp2d,gcomm_column, switchdim=.true.)
      deallocate(temp2d)
    end if

    ! column type physical state variable - zsno

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'ZSNO',['column','levsno'], &
            long_name='snow layer depth', units='m', switchdim=.true.)
    else if (flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'ZSNO') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        allocate(temp2d(begc:endc,-nlevsno+1:0))
        call clm_readvar(ncid,'ZSNO',temp2d,gcomm_column, switchdim=.true.)
        cptr%cps%z(begc:endc,-nlevsno+1:0) = temp2d(begc:endc,-nlevsno+1:0)
        deallocate(temp2d)
      end if
    else if ( flag == 'write') then
      allocate(temp2d(begc:endc,-nlevsno+1:0))
      temp2d(begc:endc,-nlevsno+1:0) = cptr%cps%z(begc:endc,-nlevsno+1:0)
      call clm_writevar(ncid,'ZSNO',temp2d,gcomm_column, switchdim=.true.)
      deallocate(temp2d)
    end if

    ! column type physical state variable - zisno

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'ZISNO',['column','levsno'], &
              long_name='snow interface depth', units='m', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'ZISNO') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        allocate(temp2d(begc:endc,-nlevsno:-1))
        call clm_readvar(ncid,'ZISNO',temp2d,gcomm_column, switchdim=.true.)
        cptr%cps%zi(begc:endc,-nlevsno:-1) = temp2d(begc:endc,-nlevsno:-1)
        deallocate(temp2d)
      end if
    else if ( flag == 'write' ) then
      allocate(temp2d(begc:endc,-nlevsno:-1))
      temp2d(begc:endc,-nlevsno:-1) = cptr%cps%zi(begc:endc,-nlevsno:-1)
      call clm_writevar(ncid,'ZISNO',temp2d,gcomm_column, switchdim=.true.)
      deallocate(temp2d)
    end if

    ! column type physical state variable - coszen

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'coszen',['column'], &
            long_name='cosine of solar zenith angle',units='unitless')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'coszen') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'coszen',cptr%cps%coszen,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'coszen',cptr%cps%coszen,gcomm_column)
    end if

    ! landunit type physical state variable - sabs_roof_dir

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'sabs_roof_dir', &
                      ['landunit','numrad  '], &
            long_name='direct solar absorbed by roof per unit '// &
                'ground area per unit incident flux',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'sabs_roof_dir') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'sabs_roof_dir',lptr%lps%sabs_roof_dir, &
                gcomm_landunit, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'sabs_roof_dir',lptr%lps%sabs_roof_dir, &
                gcomm_landunit, switchdim=.true.)
    end if

    ! landunit type physical state variable - sabs_roof_dif

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'sabs_roof_dif', &
                      ['landunit','numrad  '], &
            long_name='diffuse solar absorbed by roof per unit '// &
                'ground area per unit incident flux',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'sabs_roof_dif') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'sabs_roof_dif',lptr%lps%sabs_roof_dif, &
                gcomm_landunit, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'sabs_roof_dif',lptr%lps%sabs_roof_dif, &
                gcomm_landunit, switchdim=.true.)
    end if

    ! landunit type physical state variable - sabs_sunwall_dir

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'sabs_sunwall_dir', &
                      ['landunit','numrad  '], &
            long_name='direct solar absorbed by sunwall per unit'// &
                'wall area per unit incident flux',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'sabs_sunwall_dir') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'sabs_sunwall_dir',lptr%lps%sabs_sunwall_dir, &
                gcomm_landunit, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'sabs_sunwall_dir',lptr%lps%sabs_sunwall_dir, &
                gcomm_landunit, switchdim=.true.)
    end if

    ! landunit type physical state variable - sabs_sunwall_dif

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'sabs_sunwall_dif', &
                      ['landunit','numrad  '], &
            long_name='diffuse solar absorbed by sunwall per unit'// &
                'wall area per unit incident flux',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'sabs_sunwall_dif') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'sabs_sunwall_dif',lptr%lps%sabs_sunwall_dif, &
                gcomm_landunit, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'sabs_sunwall_dif',lptr%lps%sabs_sunwall_dif, &
                gcomm_landunit, switchdim=.true.)
    end if

    ! landunit type physical state variable - sabs_shadewall_dir

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'sabs_shadewall_dir', &
                      ['landunit','numrad  '], &
            long_name='direct solar absorbed by shadewall per unit'// &
                'wall area per unit incident flux',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'sabs_shadewall_dir') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'sabs_shadewall_dir', &
                lptr%lps%sabs_shadewall_dir, gcomm_landunit, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'sabs_shadewall_dir',lptr%lps%sabs_shadewall_dir, &
                gcomm_landunit, switchdim=.true.)
    end if

    ! landunit type physical state variable - sabs_shadewall_dif

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'sabs_shadewall_dif', &
                      ['landunit','numrad  '], &
            long_name='diffuse solar absorbed by shadewall per unit'// &
                'wall area per unit incident flux',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'sabs_shadewall_dif') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'sabs_shadewall_dif', &
                lptr%lps%sabs_shadewall_dif,gcomm_landunit, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'sabs_shadewall_dif', &
              lptr%lps%sabs_shadewall_dif,gcomm_landunit, switchdim=.true.)
    end if

    ! landunit type physical state variable - sabs_improad_dir

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'sabs_improad_dir', &
                      ['landunit','numrad  '], &
            long_name='direct solar absorbed by impervious road per unit'// &
                'ground area per unit incident flux',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'sabs_improad_dir') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'sabs_improad_dir', &
                lptr%lps%sabs_improad_dir,gcomm_landunit, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'sabs_improad_dir', &
              lptr%lps%sabs_improad_dir,gcomm_landunit, switchdim=.true.)
    end if

    ! landunit type physical state variable - sabs_improad_dif

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'sabs_improad_dif', &
                      ['landunit','numrad  '], &
            long_name='diffuse solar absorbed by impervious road per unit'// &
                'ground area per unit incident flux',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'sabs_improad_dif') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'sabs_improad_dif', &
                lptr%lps%sabs_improad_dif,gcomm_landunit, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'sabs_improad_dif', &
              lptr%lps%sabs_improad_dif,gcomm_landunit, switchdim=.true.)
    end if

    ! landunit type physical state variable - sabs_perroad_dir

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'sabs_perroad_dir', &
                      ['landunit','numrad  '], &
            long_name='direct solar absorbed by pervious road per unit '//&
                'ground area per unit incident flux',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'sabs_perroad_dir') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'sabs_perroad_dir', &
                lptr%lps%sabs_perroad_dir,gcomm_landunit, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'sabs_perroad_dir', &
              lptr%lps%sabs_perroad_dir,gcomm_landunit, switchdim=.true.)
    end if

    ! landunit type physical state variable - sabs_perroad_dif

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'sabs_perroad_dif', &
                      ['landunit','numrad  '], &
            long_name='diffuse solar absorbed by pervious road per unit '//&
                'ground area per unit incident flux',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'sabs_perroad_dif') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'sabs_perroad_dif', &
                lptr%lps%sabs_perroad_dif,gcomm_landunit, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'sabs_perroad_dif', &
              lptr%lps%sabs_perroad_dif,gcomm_landunit, switchdim=.true.)
    end if

    ! landunit type physical state variable - vf_sr

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'vf_sr', ['landunit'], &
            long_name='view factor of sky for road',units='')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'vf_sr') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'vf_sr',lptr%lps%vf_sr,gcomm_landunit)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'vf_sr',lptr%lps%vf_sr,gcomm_landunit)
    end if

    ! landunit type physical state variable - vf_wr

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'vf_wr', ['landunit'], &
            long_name='view factor of one wall for road',units='')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'vf_wr') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'vf_wr',lptr%lps%vf_wr,gcomm_landunit)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'vf_wr',lptr%lps%vf_wr,gcomm_landunit)
    end if

    ! landunit type physical state variable - vf_sw

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'vf_sw', ['landunit'], &
            long_name='view factor of sky for one wall',units='')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'vf_sw') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'vf_sw',lptr%lps%vf_sw,gcomm_landunit)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'vf_sw',lptr%lps%vf_sw,gcomm_landunit)
    end if

    ! landunit type physical state variable - vf_rw

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'vf_rw', ['landunit'], &
            long_name='view factor of road for one wall',units='')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'vf_rw') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'vf_rw',lptr%lps%vf_rw,gcomm_landunit)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'vf_rw',lptr%lps%vf_rw,gcomm_landunit)
    end if

    ! landunit type physical state variable - vf_ww

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'vf_ww', ['landunit'], &
            long_name='view factor of opposing wall for one wall',units='')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'vf_ww') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'vf_ww',lptr%lps%vf_ww,gcomm_landunit)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'vf_ww',lptr%lps%vf_ww,gcomm_landunit)
    end if

    ! landunit type physical state variable - taf

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'taf', ['landunit'], &
            long_name='urban canopy air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'taf') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'taf',lptr%lps%taf,gcomm_landunit)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'taf',lptr%lps%taf,gcomm_landunit)
    end if

    ! landunit type physical state variable - qaf

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'qaf', ['landunit'], &
            long_name='urban canopy specific humidity',units='kg/kg')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'qaf') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'qaf',lptr%lps%qaf,gcomm_landunit)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'qaf',lptr%lps%qaf,gcomm_landunit)
    end if

    ! pft type physical state variable - albd

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albd', ['pft   ','numrad'], &
            long_name='surface albedo (direct) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'albd') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        if ( rcmtimer%start( ) ) then
          do_initsurfalb = .true.
        else
          call clm_readvar(ncid,'albd',pptr%pps%albd,gcomm_pft, switchdim=.true.)
        end if
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albd',pptr%pps%albd,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - albi

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albi', ['pft   ','numrad'], &
            long_name='surface albedo (diffuse) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'albi') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        if ( rcmtimer%start( ) ) then
          do_initsurfalb = .true.
        else
          call clm_readvar(ncid,'albi',pptr%pps%albi,gcomm_pft, switchdim=.true.)
        end if
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albi',pptr%pps%albi,gcomm_pft, switchdim=.true.)
    end if

    ! column type physical state variable - albgrd

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albgrd', ['column','numrad'], &
            long_name='ground albedo (direct) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'albgrd') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'albgrd',cptr%cps%albgrd,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albgrd',cptr%cps%albgrd,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - albgri

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albgri', ['column','numrad'], &
            long_name='ground albedo (diffuse) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'albgri') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        call clm_readvar(ncid,'albgri',cptr%cps%albgri,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albgri',cptr%cps%albgri,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - albsod

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albsod', ['column','numrad'], &
            long_name='soil albedo (direct) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'albsod') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        if ( rcmtimer%start( ) ) then
          do_initsurfalb = .true.
        else
          call clm_readvar(ncid,'albsod',cptr%cps%albsod,gcomm_column, switchdim=.true.)
        end if
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albsod',cptr%cps%albsod,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - albsoi

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albsoi', ['column','numrad'], &
            long_name='soil albedo (diffuse) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'albsoi') ) then
        call fatal(__FILE__,__LINE__,'clm_now_stopping')
      else
        if ( rcmtimer%start( ) ) then
          do_initsurfalb = .true.
        else
          call clm_readvar(ncid,'albsoi',cptr%cps%albsoi,gcomm_column, switchdim=.true.)
        end if
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albsoi',cptr%cps%albsoi,gcomm_column, switchdim=.true.)
    end if

#ifdef SNICAR_FRC
    ! column type physical state variable - albgrd_bc

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albgrd_bc', ['column','numrad'], &
            long_name='ground albedo without BC (direct) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'albgrd_bc') ) then
        if (myid == italk) then
          write(stderr,*) &
            "SNICAR: can't find albgrd_bc in restart (or initial) file..."
          write(stderr,*) "Initialize albgrd_bc to albgrd"
        end if
        do c = begc , endc
          cptr%cps%albgrd_bc(c,:) = cptr%cps%albgrd(c,:)
        end do
      else
        call clm_readvar(ncid,'albgrd_bc',cptr%cps%albgrd_bc,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albgrd_bc',cptr%cps%albgrd_bc,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - albgri_bc

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albgri_bc', ['column','numrad'], &
            long_name='ground albedo without BC (diffuse) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'albgri_bc') ) then
        if (myid == italk) then
          write(stderr,*) &
            "SNICAR: can't find albgri_bc in restart (or initial) file..."
          write(stderr,*) "Initialize albgri_bc to albgri"
        end if
        do c = begc , endc
          cptr%cps%albgri_bc(c,:) = cptr%cps%albgri(c,:)
        end do
      else
        call clm_readvar(ncid,'albgri_bc',cptr%cps%albgri_bc,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albgri_bc',cptr%cps%albgri_bc,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - albgrd_pur

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albgrd_pur',['column','numrad'], &
            long_name='pure snow ground albedo (direct) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'albgrd_pur') ) then
        if (myid == italk) then
          write(stderr,*) &
            "SNICAR: can't find albgrd_pur in restart (or initial) file..."
          write(stderr,*) "Initialize albgrd_pur to albgrd"
        end if
        do c = begc , endc
          cptr%cps%albgrd_pur(c,:) = cptr%cps%albgrd(c,:)
        end do
      else
        call clm_readvar(ncid,'albgrd_pur',cptr%cps%albgrd_pur,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albgrd_pur',cptr%cps%albgrd_pur,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - albgri_pur

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albgri_pur',['column','numrad'], &
            long_name='pure snow ground albedo (diffuse) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'albgri_pur') ) then
        if (myid == italk) then
          write(stderr,*) &
            "SNICAR: can't find albgri_pur in restart (or initial) file..."
          write(stderr,*) "Initialize albgri_pur to albgrd"
        end if
        do c = begc , endc
          cptr%cps%albgri_pur(c,:) = cptr%cps%albgri(c,:)
        end do
      else
        call clm_readvar(ncid,'albgri_pur',cptr%cps%albgri_pur,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albgri_pur',cptr%cps%albgri_pur,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - albgrd_oc

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albgrd_oc', ['column','numrad'], &
            long_name='ground albedo without OC (direct) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'albgrd_oc') ) then
        if (myid == italk) then
          write(stderr,*) &
            "SNICAR: can't find albgrd_oc in restart (or initial) file..."
          write(stderr,*) "Initialize albgrd_oc to albgrd"
        end if
        do c = begc , endc
          cptr%cps%albgrd_oc(c,:) = cptr%cps%albgrd(c,:)
        end do
      else
        call clm_readvar(ncid,'albgrd_oc',cptr%cps%albgrd_oc,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albgrd_oc',cptr%cps%albgrd_oc,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - albgri_oc

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albgri_oc', ['column','numrad'], &
            long_name='ground albedo without OC (diffuse) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'albgri_oc') ) then
        if (myid == italk) then
          write(stderr,*) &
            "SNICAR: can't find albgri_oc in restart (or initial) file..."
          write(stderr,*) "Initialize albgri_oc to albgri"
        end if
        do c = begc , endc
          cptr%cps%albgri_oc(c,:) = cptr%cps%albgri(c,:)
        end do
      else
        call clm_readvar(ncid,'albgri_oc',cptr%cps%albgri_oc,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albgri_oc',cptr%cps%albgri_oc,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - albgrd_dst

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albgrd_dst', ['column','numrad'], &
            long_name='ground albedo without dust (direct) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'albgrd_dst') ) then
        if (myid == italk) then
          write(stderr,*) &
            "SNICAR: can't find albgrd_dst in restart (or initial) file..."
          write(stderr,*) "Initialize albgrd_dst to albgrd"
        end if
        do c = begc , endc
          cptr%cps%albgrd_dst(c,:) = cptr%cps%albgrd(c,:)
        end do
      else
        call clm_readvar(ncid,'albgrd_dst',cptr%cps%albgrd_dst,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albgrd_dst',cptr%cps%albgrd_dst,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - albgri_dst

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albgri_dst', ['column','numrad'], &
           long_name='ground albedo without dust (diffuse) (0 to 1)',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'albgri_dst') ) then
        if (myid == italk) then
          write(stderr,*) &
            "SNICAR: can't find albgri_dst in restart (or initial) file..."
          write(stderr,*) "Initialize albgri_dst to albgrd"
        end if
        do c = begc , endc
          cptr%cps%albgri_dst(c,:) = cptr%cps%albgri(c,:)
        end do
      else
        call clm_readvar(ncid,'albgri_dst',cptr%cps%albgri_dst,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albgri_dst',cptr%cps%albgri_dst,gcomm_column, switchdim=.true.)
    end if
#endif

    ! column water state variable - h2osfc

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'H2OSFC', ['column'], &
            long_name='surface water',units='kg/m2')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'H2OSFC') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          cptr%cws%h2osfc(begc:endc) = 0.0_rk8
        end if
      else
        call clm_readvar(ncid,'H2OSFC',cptr%cws%h2osfc,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'H2OSFC',cptr%cws%h2osfc,gcomm_column)
    end if

    ! column type physical state variable - frac_h2osfc

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'FH2OSFC', ['column'], &
            long_name='fraction of ground covered by h2osfc (0 to 1)',units='')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'FH2OSFC') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          cptr%cps%frac_h2osfc(begc:endc) = 0.0_rk8
        end if
      else
        call clm_readvar(ncid,'FH2OSFC',cptr%cps%frac_h2osfc,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'FH2OSFC',cptr%cps%frac_h2osfc,gcomm_column)
    end if

   ! column energy state variable - t_h2osfc

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'TH2OSFC', ['column'], &
            long_name='surface water temperature',units='K')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'TH2OSFC') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          cptr%ces%t_h2osfc(begc:endc) = 274.0_rk8
        end if
      else
        call clm_readvar(ncid,'TH2OSFC',cptr%ces%t_h2osfc,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'TH2OSFC',cptr%ces%t_h2osfc,gcomm_column)
    end if

   ! column water state variable - h2osno

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'H2OSNO', ['column'], &
            long_name='snow water',units='kg/m2')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'H2OSNO') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'H2OSNO',cptr%cws%h2osno,gcomm_column)
      end if
      where ( cptr%cws%h2osno < 1.0e-40_rk8 )
        cptr%cws%h2osno = 0.0_rk8
      end where
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'H2OSNO',cptr%cws%h2osno,gcomm_column)
    end if

    ! column water state variable - h2osoi_liq

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'H2OSOI_LIQ', ['column','levtot'], &
            long_name='liquid water',units='kg/m2', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'H2OSOI_LIQ') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'H2OSOI_LIQ',cptr%cws%h2osoi_liq,gcomm_column, switchdim=.true.)
        where ( cptr%cws%h2osoi_liq < 1.0e-40_rk8 )
          cptr%cws%h2osoi_liq = 0.0_rk8
        end where
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'H2OSOI_LIQ',cptr%cws%h2osoi_liq,gcomm_column, switchdim=.true.)
    end if

    ! column water state variable - h2osoi_ice

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'H2OSOI_ICE', ['column','levtot'], &
            long_name='ice lens',units='kg/m2', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'H2OSOI_ICE') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'H2OSOI_ICE',cptr%cws%h2osoi_ice,gcomm_column, switchdim=.true.)
        where ( cptr%cws%h2osoi_ice < 1.0e-40_rk8 )
          cptr%cws%h2osoi_ice = 0.0_rk8
        end where
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'H2OSOI_ICE',cptr%cws%h2osoi_ice,gcomm_column, switchdim=.true.)
    end if

   ! column energy state variable - t_grnd

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_GRND', ['column'], &
            long_name='ground temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_GRND') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_GRND',cptr%ces%t_grnd,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_GRND',cptr%ces%t_grnd,gcomm_column)
    end if

   ! column urban energy state variable - eflx_urban_ac

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'URBAN_AC', ['column'], &
            long_name='urban air conditioning flux',units='W/m^2')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'URBAN_AC') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'URBAN_AC',cptr%cef%eflx_urban_ac,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'URBAN_AC',cptr%cef%eflx_urban_ac,gcomm_column)
    end if

   ! column urban energy state variable - eflx_urban_heat

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'URBAN_HEAT', ['column'], &
            long_name='urban heating flux',units='W/m^2')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'URBAN_HEAT') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'URBAN_HEAT', &
                cptr%cef%eflx_urban_heat,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'URBAN_HEAT',cptr%cef%eflx_urban_heat,gcomm_column)
    end if

   ! pft energy state variable - t_ref2m_min

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_MIN', ['pft'], &
            long_name='daily minimum of average 2 m height '// &
               'surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_MIN') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_MIN',pptr%pes%t_ref2m_min,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_MIN',pptr%pes%t_ref2m_min,gcomm_pft)
    end if

   ! pft energy state variable - t_ref2m_max

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_MAX', ['pft'], &
            long_name='daily maximum of average 2 m height '// &
               'surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_MAX') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_MAX',pptr%pes%t_ref2m_max,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_MAX',pptr%pes%t_ref2m_max,gcomm_pft)
    end if

   ! pft energy state variable - t_ref2m_min_inst

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_MIN_INST', ['pft'], &
            long_name='instantaneous daily minimum of average 2 m height '// &
               'surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_MIN_INST') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_MIN_INST', &
                pptr%pes%t_ref2m_min_inst,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_MIN_INST', &
              pptr%pes%t_ref2m_min_inst,gcomm_pft)
    end if

   ! pft energy state variable - t_ref2m_max_inst

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_MAX_INST', ['pft'], &
            long_name='instantaneous daily maximum of average 2 m height '// &
               'surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_MAX_INST') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_MAX_INST', &
                pptr%pes%t_ref2m_max_inst,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_MAX_INST', &
              pptr%pes%t_ref2m_max_inst,gcomm_pft)
    end if

    ! pft energy state variable - t_ref2m_u

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_U', ['pft'], &
            long_name='Urban 2m height surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_U') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_U',pptr%pes%t_ref2m_u,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_U',pptr%pes%t_ref2m_u,gcomm_pft)
    end if

   ! column energy state variable - t_grnd_u

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_GRND_U', ['column'], &
            long_name='Urban ground temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_GRND_U') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_GRND_U',cptr%ces%t_grnd_u,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_GRND_U',cptr%ces%t_grnd_u,gcomm_column)
    end if

   ! pft energy state variable - t_ref2m_min_u

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_MIN_U', ['pft'], &
            long_name='urban daily minimum of average 2 m height '// &
               'surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_MIN_U') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_MIN_U',pptr%pes%t_ref2m_min_u,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_MIN_U',pptr%pes%t_ref2m_min_u,gcomm_pft)
    end if

   ! pft energy state variable - t_ref2m_max_u

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_MAX_U', ['pft'], &
            long_name='urban daily maximum of average 2 m height '// &
               'surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_MAX_U') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_MAX_U',pptr%pes%t_ref2m_max_u,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_MAX_U',pptr%pes%t_ref2m_max_u,gcomm_pft)
    end if

   ! pft energy state variable - t_ref2m_min_inst_u

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_MIN_INST_U', ['pft'], &
            long_name='urban instantaneous daily min of average 2 m height '// &
               'surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_MIN_INST_U') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_MIN_INST_U', &
                pptr%pes%t_ref2m_min_inst_u,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_MIN_INST_U', &
              pptr%pes%t_ref2m_min_inst_u,gcomm_pft)
    end if

   ! pft energy state variable - t_ref2m_max_inst_u

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_MAX_INST_U', ['pft'], &
            long_name='urban instantaneous daily max of average 2 m height '// &
               'surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_MAX_INST_U') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_MAX_INST_U', &
                pptr%pes%t_ref2m_max_inst_u,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_MAX_INST_U', &
              pptr%pes%t_ref2m_max_inst_u,gcomm_pft)
    end if

    ! pft energy state variable - t_ref2m_r

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_R', ['pft'], &
            long_name='Rural 2m height surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_R') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_R',pptr%pes%t_ref2m_r,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_R',pptr%pes%t_ref2m_r,gcomm_pft)
    end if

   ! column energy state variable - t_grnd_r

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_GRND_R', ['column'], &
            long_name='Rural ground temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_GRND_R') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_GRND_R',cptr%ces%t_grnd_r,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_GRND_R',cptr%ces%t_grnd_r,gcomm_column)
    end if

   ! pft energy state variable - t_ref2m_min_r

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_MIN_R', ['pft'], &
            long_name='rural daily minimum of average 2 m height '// &
            'surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_MIN_R') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_MIN_R',pptr%pes%t_ref2m_min_r,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_MIN_R',pptr%pes%t_ref2m_min_r,gcomm_pft)
    end if

   ! pft energy state variable - t_ref2m_max_r

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_MAX_R', ['pft'], &
            long_name='rural daily maximum of average 2 m height '// &
            'surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_MAX_R') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_MAX_R',pptr%pes%t_ref2m_max_r,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_MAX_R',pptr%pes%t_ref2m_max_r,gcomm_pft)
    end if

   ! pft energy state variable - t_ref2m_min_inst_r

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_MIN_INST_R', ['pft'], &
            long_name='rural instantaneous daily min of average 2 m height '// &
            'surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_MIN_INST_R') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_MIN_INST_R', &
                pptr%pes%t_ref2m_min_inst_r,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_MIN_INST_R', &
              pptr%pes%t_ref2m_min_inst_r,gcomm_pft)
    end if

   ! pft energy state variable - t_ref2m_max_inst_r

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M_MAX_INST_R', ['pft'], &
            long_name='rural instantaneous daily max of average 2 m height '// &
            'surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_REF2M_MAX_INST_R') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_REF2M_MAX_INST_R', &
                pptr%pes%t_ref2m_max_inst_r,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M_MAX_INST_R', &
              pptr%pes%t_ref2m_max_inst_r,gcomm_pft)
    end if

    ! column energy state variable - t_soisno

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_SOISNO', ['column','levtot'], &
            long_name='soil-snow temperature',units='K', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_SOISNO') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_SOISNO',cptr%ces%t_soisno,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_SOISNO',cptr%ces%t_soisno,gcomm_column, switchdim=.true.)
    end if

    ! column type energy state variable - t_lake

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_LAKE', ['column','levlak'], &
            long_name='lake temperature',units='K', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_LAKE') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_LAKE',cptr%ces%t_lake,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_LAKE',cptr%ces%t_lake,gcomm_column, switchdim=.true.)
    end if

    ! pft physical state variable - frac_veg_nosno_alb

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_integer,ncid,'FRAC_VEG_NOSNO_ALB', ['pft'], &
            long_name='fraction of vegetation not covered by snow',units='')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'FRAC_VEG_NOSNO_ALB') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'FRAC_VEG_NOSNO_ALB', &
                pptr%pps%frac_veg_nosno_alb,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'FRAC_VEG_NOSNO_ALB', &
              pptr%pps%frac_veg_nosno_alb,gcomm_pft)
    end if

    ! pft type physical state variable - fwet

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'FWET', ['pft'], &
            long_name='fraction of canopy that is wet',units='')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'FWET') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'FWET',pptr%pps%fwet,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'FWET',pptr%pps%fwet,gcomm_pft)
    end if

    ! pft type physical state variable - tlai

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'tlai', ['pft'], &
            long_name='one-sided leaf area index, no burying by snow',units='')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'tlai') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'tlai',pptr%pps%tlai,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'tlai',pptr%pps%tlai,gcomm_pft)
    end if

    ! pft type physical state variable - tsai

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'tsai', ['pft'], &
            long_name='one-sided stem area index, no burying by snow',units='')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'tsai') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'tsai',pptr%pps%tsai,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'tsai',pptr%pps%tsai,gcomm_pft)
    end if

    ! pft type physical state variable - tlai_z

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'tlai_z', ['pft   ','levcan'], &
            long_name='tlai increment for canopy layer',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'tsai') ) then
        if ( myid == italk ) then
          write(stderr,*) "can't find tlai_z in restart (or initial) file..."
          write(stderr,*) "Initialize tlai_z to tlai/nlevcan"
        end if
        do p = begp , endp
          do iv = 1 , nlevcan
            pptr%pps%tlai_z(p,iv) = pptr%pps%tlai(p)/nlevcan
          end do
        end do
      else
        call clm_readvar(ncid,'tlai_z',pptr%pps%tlai_z,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'tlai_z',pptr%pps%tlai_z,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - tsai_z

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'tsai_z', ['pft   ','levcan'], &
            long_name='tsai increment for canopy layer',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'tsai') ) then
        if ( myid == italk ) then
          write(stderr,*) "can't find tsai_z in restart (or initial) file..."
          write(stderr,*) "Initialize tsai_z to tsai/nlevcan"
        end if
        do p = begp , endp
          do iv = 1 , nlevcan
            pptr%pps%tsai_z(p,iv) = pptr%pps%tsai(p)/nlevcan
          end do
        end do
      else
        call clm_readvar(ncid,'tsai_z',pptr%pps%tsai_z,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'tsai_z',pptr%pps%tsai_z,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - ncan

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_integer,ncid,'ncan', ['pft'], &
            long_name='number of canopy layer',units='')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'ncan') ) then
        if ( myid == italk ) then
          write(stderr,*) "can't find ncan in restart (or initial) file..."
          write(stderr,*) "Initialize ncan to nlevcan"
        end if
        do p = begp , endp
          pptr%pps%ncan(p) = nlevcan
        end do
      else
        call clm_readvar(ncid,'ncan',pptr%pps%ncan,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'ncan',pptr%pps%ncan,gcomm_pft)
    end if

    ! pft type physical state variable - nrad

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_integer,ncid,'nrad', ['pft'], &
            long_name='number of canopy layer above snow',units='')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'nrad') ) then
        if ( myid == italk ) then
          write(stderr,*) "can't find nrad in restart (or initial) file..."
          write(stderr,*) "Initialize nrad to nlevcan"
        end if
        do p = begp , endp
          pptr%pps%nrad(p) = nlevcan
        end do
      else
        call clm_readvar(ncid,'nrad',pptr%pps%nrad,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'nrad',pptr%pps%nrad,gcomm_pft)
    end if

    ! pft type physical state variable - mlaidiff

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'mlaidiff', ['pft'], &
            long_name='difference between lai month one and month two',units='')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'mlaidiff') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'mlaidiff',pptr%pps%mlaidiff,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'mlaidiff',pptr%pps%mlaidiff,gcomm_pft)
    end if

    ! pft type physical state variable - elai

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'elai', ['pft'], &
            long_name='one-sided leaf area index, with burying by snow', &
            units='')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'elai') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'elai',pptr%pps%elai,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'elai',pptr%pps%elai,gcomm_pft)
    end if

    ! pft type physical state variable - esai

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'esai', ['pft'], &
            long_name='one-sided stem area index, with burying by snow', &
            units='')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'esai') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'esai',pptr%pps%esai,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'esai',pptr%pps%esai,gcomm_pft)
    end if

    ! pft type physical state variable - fsun

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fsun', ['pft'], &
            long_name='sunlit fraction of canopy',units='')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'fsun') ) then
        if ( rcmtimer%start( ) ) then
          do p = begp , endp
            if ( is_nan( pptr%pps%fsun(p) ) )then
              pptr%pps%fsun(p) = spval
            end if
          end do
        else
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'fsun',pptr%pps%fsun,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fsun',pptr%pps%fsun,gcomm_pft)
    end if

    ! pft type physical state variable - vcmaxcintsun

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'vcmaxcintsun', ['pft'], &
            long_name='sunlit canopy scaling coefficient',units='')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'vcmaxcintsun') ) then
        write(stderr,*) "can't find vcmaxcintsun in restart (or init) file..."
        write(stderr,*) "Initialize vcmaxcintsun to 1"
        do p = begp , endp
          pptr%pps%vcmaxcintsun(p) = 1._rk8
        end do
      else
        call clm_readvar(ncid,'vcmaxcintsun',pptr%pps%vcmaxcintsun,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'vcmaxcintsun',pptr%pps%vcmaxcintsun,gcomm_pft)
    end if

    ! pft type physical state variable - vcmaxcintsha

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'vcmaxcintsha', ['pft'], &
            long_name='shaded canopy scaling coefficient',units='')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'vcmaxcintsha') ) then
        write(stderr,*) "can't find vcmaxcintsha in restart (or init) file..."
        write(stderr,*) "Initialize vcmaxcintsha to 1"
        do p = begp , endp
          pptr%pps%vcmaxcintsha(p) = 1._rk8
        end do
      else
        call clm_readvar(ncid,'vcmaxcintsha',pptr%pps%vcmaxcintsha,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'vcmaxcintsha',pptr%pps%vcmaxcintsha,gcomm_pft)
    end if

    ! pft type physical state variable - htop

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'htop', ['pft'], &
            long_name='canopy top',units='m')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'htop') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'htop',pptr%pps%htop,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'htop',pptr%pps%htop,gcomm_pft)
    end if

    ! pft type physical state variable - hbot

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'hbot', ['pft'], &
            long_name='canopy botton',units='m')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'hbot') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'hbot',pptr%pps%hbot,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'hbot',pptr%pps%hbot,gcomm_pft)
    end if

    ! pft type physical state variable - fabd

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fabd', ['pft   ','numrad'], &
            long_name='flux absorbed by veg per unit direct flux',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'fabd') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'fabd',pptr%pps%fabd,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fabd',pptr%pps%fabd,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - fabi

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fabi', ['pft   ','numrad'], &
            long_name='flux absorbed by veg per unit diffuse flux',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'fabi') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'fabi',pptr%pps%fabi,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fabi',pptr%pps%fabi,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - fabd_sun

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fabd_sun', ['pft   ','numrad'], &
            long_name='flux absorbed by sunlit leaf per unit direct flux', &
            units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'fabd_sun') ) then
        write(stderr,*) "can't find fabd_sun in restart (or initial) file..."
        write(stderr,*) "Initialize fabd_sun to fabd/2"
        do p = begp , endp
          pptr%pps%fabd_sun(p,:) = pptr%pps%fabd(p,:)/2._rk8
        end do
      else
        call clm_readvar(ncid,'fabd_sun',pptr%pps%fabd_sun,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fabd_sun',pptr%pps%fabd_sun,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - fabd_sha

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fabd_sha', ['pft   ','numrad'], &
            long_name='flux absorbed by shaded leaf per unit direct flux', &
            units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'fabd_sha') ) then
        write(stderr,*) "can't find fabd_sha in restart (or initial) file..."
        write(stderr,*) "Initialize fabd_sha to fabd/2"
        do p = begp , endp
          pptr%pps%fabd_sha(p,:) = pptr%pps%fabd(p,:)/2._rk8
        end do
      else
        call clm_readvar(ncid,'fabd_sha',pptr%pps%fabd_sha,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fabd_sha',pptr%pps%fabd_sha,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - fabi_sun

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fabi_sun', ['pft   ','numrad'], &
            long_name='flux absorbed by sunlit leaf per unit diffuse flux', &
            units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'fabi_sun') ) then
        write(stderr,*) "can't find fabi_sun in restart (or initial) file..."
        write(stderr,*) "Initialize fabi_sun to fabi/2"
        do p = begp , endp
          pptr%pps%fabi_sun(p,:) = pptr%pps%fabi(p,:)/2._rk8
        end do
      else
        call clm_readvar(ncid,'fabi_sun',pptr%pps%fabi_sun,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fabi_sun',pptr%pps%fabi_sun,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - fabi_sha

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fabi_sha', ['pft   ','numrad'], &
            long_name='flux absorbed by shaded leaf per unit diffuse flux', &
            units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'fabi_sha') ) then
        write(stderr,*) "can't find fabi_sha in restart (or initial) file..."
        write(stderr,*) "Initialize fabi_sha to fabi/2"
        do p = begp , endp
          pptr%pps%fabi_sha(p,:) = pptr%pps%fabi(p,:)/2._rk8
        end do
      else
        call clm_readvar(ncid,'fabi_sha',pptr%pps%fabi_sha,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fabi_sha',pptr%pps%fabi_sha,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - fabd_sun_z

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fabd_sun_z', ['pft   ','levcan'], &
            long_name='absorbed sunlit leaf direct PAR (per unit lai+sai) '// &
                  'for canopy layer',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'fabd_sun_z') ) then
        write(stderr,*) "can't find fabd_sun_z in restart (or initial) file..."
        write(stderr,*) "Initialize fabd_sun_z to fabd/2/nlevcan"
        do p = begp , endp
          do iv = 1 , nlevcan
            pptr%pps%fabd_sun_z(p,iv) = (pptr%pps%fabd(p,1)/2._rk8)/real(nlevcan,rk8)
          end do
        end do
      else
        call clm_readvar(ncid,'fabd_sun_z',pptr%pps%fabd_sun_z,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fabd_sun_z',pptr%pps%fabd_sun_z,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - fabd_sha_z

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fabd_sha_z', ['pft   ','levcan'], &
            long_name='absorbed shaded leaf direct PAR (per unit lai+sai) '// &
                  'for canopy layer',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'fabd_sha_z') ) then
        write(stderr,*) "can't find fabd_sha_z in restart (or initial) file..."
        write(stderr,*) "Initialize fabd_sha_z to fabd/2/nlevcan"
        do p = begp , endp
          do iv = 1 , nlevcan
            pptr%pps%fabd_sha_z(p,iv) = (pptr%pps%fabd(p,1)/2._rk8)/real(nlevcan,rk8)
          end do
        end do
      else
        call clm_readvar(ncid,'fabd_sha_z',pptr%pps%fabd_sha_z,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fabd_sha_z',pptr%pps%fabd_sha_z,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - fabi_sun_z

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fabi_sun_z', ['pft   ','levcan'], &
            long_name='absorbed sunlit leaf diffuse PAR (per unit lai+sai) '// &
                  'for canopy layer',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'fabi_sun_z') ) then
        write(stderr,*) "can't find fabi_sun_z in restart (or initial) file..."
        write(stderr,*) "Initialize fabi_sun_z to fabi/2/nlevcan"
        do p = begp , endp
          do iv = 1 , nlevcan
            pptr%pps%fabi_sun_z(p,iv) = (pptr%pps%fabi(p,1)/2._rk8)/real(nlevcan,rk8)
          end do
        end do
      else
        call clm_readvar(ncid,'fabi_sun_z',pptr%pps%fabi_sun_z,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fabi_sun_z',pptr%pps%fabi_sun_z,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - fabi_sha_z

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fabi_sha_z', ['pft   ','levcan'], &
            long_name='absorbed shaded leaf diffuse PAR (per unit lai+sai) '// &
                  'for canopy layer',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'fabi_sha_z') ) then
        write(stderr,*) "can't find fabi_sha_z in restart (or initial) file..."
        write(stderr,*) "Initialize fabi_sha_z to fabi/2/nlevcan"
        do p = begp , endp
          do iv = 1 , nlevcan
            pptr%pps%fabi_sha_z(p,iv) = (pptr%pps%fabi(p,1)/2._rk8)/real(nlevcan,rk8)
          end do
        end do
      else
        call clm_readvar(ncid,'fabi_sha_z',pptr%pps%fabi_sha_z,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fabi_sha_z',pptr%pps%fabi_sha_z,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - fsun_z

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'fsun_z', ['pft   ','levcan'], &
            long_name='absorbed shaded leaf diffuse PAR (per unit lai+sai) '// &
                  'for canopy layer',units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'fsun_z') ) then
        write(stderr,*) "can't find fsun_z in restart (or initial) file..."
        write(stderr,*) "Initialize fsun_z to 0"
        do p = begp , endp
          do iv = 1 , nlevcan
            pptr%pps%fsun_z(p,iv) = 0._rk8
          end do
        end do
      else
        call clm_readvar(ncid,'fsun_z',pptr%pps%fsun_z,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'fsun_z',pptr%pps%fsun_z,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - ftdd

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'ftdd', ['pft   ','numrad'], &
              long_name='down direct flux below veg per unit direct flux', &
              units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'ftdd') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'ftdd',pptr%pps%ftdd,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'ftdd',pptr%pps%ftdd,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - ftid

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'ftid', ['pft   ','numrad'], &
              long_name='down diffuse flux below veg per unit direct flux', &
              units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'ftid') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'ftid',pptr%pps%ftid,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'ftid',pptr%pps%ftid,gcomm_pft, switchdim=.true.)
    end if

    ! pft type physical state variable - ftii

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'ftii', ['pft   ','numrad'], &
              long_name='down diffuse flux below veg per unit diffuse flux', &
              units='', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'ftii') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'ftii',pptr%pps%ftii,gcomm_pft, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'ftii',pptr%pps%ftii,gcomm_pft, switchdim=.true.)
    end if

    ! pft energy state variable - t_veg

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_VEG', ['pft'], &
              long_name='vegetation temperature',units='K')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'T_VEG') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_VEG',pptr%pes%t_veg,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_VEG',pptr%pes%t_veg,gcomm_pft)
    end if

    ! pft energy state variable - t_ref2m

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'T_REF2M', ['pft'], &
              long_name='2m height surface air temperature',units='K')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'T_REF2M') ) then
        if ( rcmtimer%start( ) ) then
          if ( allocate_all_vegpfts ) then
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'T_REF2M',pptr%pes%t_ref2m,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'T_REF2M',pptr%pes%t_ref2m,gcomm_pft)
    end if

    ! pft type water state variable - h2ocan

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'H2OCAN', ['pft'], &
              long_name='canopy water',units='kg/m2')
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'H2OCAN') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'H2OCAN',pptr%pws%h2ocan,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'H2OCAN',pptr%pws%h2ocan,gcomm_pft)
    end if

   ! column irrigation variable - n_irrig_steps_left

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_integer,ncid,'n_irrig_steps_left', ['column'], &
              long_name='number of irrigation time steps left',units='#')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'n_irrig_steps_left') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          cptr%cps%n_irrig_steps_left = 0
        end if
      else
        call clm_readvar(ncid,'n_irrig_steps_left', &
                cptr%cps%n_irrig_steps_left,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'n_irrig_steps_left', &
              cptr%cps%n_irrig_steps_left,gcomm_column)
    end if

   ! column irrigation variable - irrig_rate

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'irrig_rate', ['column'], &
              long_name='irrigation rate',units='mm/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'irrig_rate') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          cptr%cps%irrig_rate = 0.0_rk8
        end if
      else
        call clm_readvar(ncid,'irrig_rate',cptr%cps%irrig_rate,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'irrig_rate',cptr%cps%irrig_rate,gcomm_column)
    end if

    ! ------------------------------------------------------------
    ! Determine volumetric soil water (for read only)
    ! ------------------------------------------------------------

    if ( flag == 'read' ) then
      do c = begc , endc
        l = clandunit(c)
        if ( ctype(c) == icol_sunwall .or. &
             ctype(c) == icol_shadewall .or. &
             ctype(c) == icol_roof ) then
          nlevs = nlevurb
        else
          nlevs = nlevgrnd
        end if
        ! NOTE: THIS IS A MEMORY INEFFICIENT COPY
        if ( ltype(l) /= istdlak ) then
          ! This calculation is now done for lakes in initSLake.
          do j = 1 , nlevs
            if ( cptr%cws%h2osoi_liq(c,j) < 1.0e-40_rk8 ) then
              cptr%cws%h2osoi_liq(c,j) = 0.0_rk8
            end if
            if ( cptr%cws%h2osoi_ice(c,j) < 1.0e-40_rk8 ) then
              cptr%cws%h2osoi_ice(c,j) = 0.0_rk8
            end if
            cptr%cws%h2osoi_vol(c,j) = &
                    cptr%cws%h2osoi_liq(c,j)/(cptr%cps%dz(c,j)*denh2o) + &
                    cptr%cws%h2osoi_ice(c,j)/(cptr%cps%dz(c,j)*denice)
          end do
        end if
      end do

      ! ------------------------------------------------------------
      ! If initial run -- ensure that water is properly bounded
      ! ------------------------------------------------------------

      if ( rcmtimer%start( ) ) then
        do c = begc , endc
          l = clandunit(c)
          if ( ctype(c) == icol_sunwall .or. &
               ctype(c) == icol_shadewall .or. &
               ctype(c) == icol_roof ) then
            nlevs = nlevurb
          else
            nlevs = nlevgrnd
          end if
          do j = 1 , nlevs
            l = clandunit(c)
            if ( ltype(l) == istsoil .or. ltype(l) == istcrop ) then
              cptr%cws%h2osoi_liq(c,j) = max(0._rk8,cptr%cws%h2osoi_liq(c,j))
              cptr%cws%h2osoi_ice(c,j) = max(0._rk8,cptr%cws%h2osoi_ice(c,j))
              cptr%cws%h2osoi_vol(c,j) = &
                      cptr%cws%h2osoi_liq(c,j)/(cptr%cps%dz(c,j)*denh2o) + &
                      cptr%cws%h2osoi_ice(c,j)/(cptr%cps%dz(c,j)*denice)
              if ( j == 1 ) then
                maxwatsat = (cptr%cps%watsat(c,j)*cptr%cps%dz(c,j)*1000.0_rk8 + &
                        pondmx)/(cptr%cps%dz(c,j)*1000.0_rk8)
              else
                maxwatsat = cptr%cps%watsat(c,j)
              end if
              if ( cptr%cws%h2osoi_vol(c,j) > maxwatsat ) then
                excess = (cptr%cws%h2osoi_vol(c,j) - maxwatsat) * &
                        cptr%cps%dz(c,j)*1000.0_rk8
                totwat = cptr%cws%h2osoi_liq(c,j) + cptr%cws%h2osoi_ice(c,j)
                cptr%cws%h2osoi_liq(c,j) = cptr%cws%h2osoi_liq(c,j) - &
                        (cptr%cws%h2osoi_liq(c,j)/totwat) * excess
                cptr%cws%h2osoi_ice(c,j) = cptr%cws%h2osoi_ice(c,j) - &
                        (cptr%cws%h2osoi_ice(c,j)/totwat) * excess
              end if
              cptr%cws%h2osoi_liq(c,j) = max(watmin,cptr%cws%h2osoi_liq(c,j))
              cptr%cws%h2osoi_ice(c,j) = max(watmin,cptr%cws%h2osoi_ice(c,j))
              cptr%cws%h2osoi_vol(c,j) = &
                      cptr%cws%h2osoi_liq(c,j)/(cptr%cps%dz(c,j)*denh2o) + &
                      cptr%cws%h2osoi_ice(c,j)/(cptr%cps%dz(c,j)*denice)
            end if
          end do
        end do
      end if
    end if
    !
    ! Perturb initial conditions if not restart
    !
    if ( .not. lstart .and. flag=='read' .and. pertlim /= 0.0_rk8 ) then
      call perturbIC( lptr )
    end if
    !
    ! variables needed for SNICAR
    !

    ! column type physical state variable - snw_rds

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'snw_rds', ['column','levsno'], &
              long_name='snow layer effective radius',units='um', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'snw_rds') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          write(stderr,*) 'SNICAR: This is an initial run (not a restart)'
          write(stderr,*) 'grain size/aerosol mass data are not defined in &
                         &initial condition file'
          write(stderr,*) 'Initialize snow effective radius to fresh snow &
                         &value , and snow/aerosol masses to zero.'
          do c = begc , endc
            if ( cptr%cps%snl(c) < 0 ) then
              cptr%cps%snw_rds(c,cptr%cps%snl(c)+1:0) = snw_rds_min
              cptr%cps%snw_rds(c,-nlevsno+1:cptr%cps%snl(c)) = 0._rk8
              cptr%cps%snw_rds_top(c) = snw_rds_min
              cptr%cps%sno_liq_top(c) = &
                      cptr%cws%h2osoi_liq(c,cptr%cps%snl(c)+1) / &
                      (cptr%cws%h2osoi_liq(c,cptr%cps%snl(c)+1) + &
                       cptr%cws%h2osoi_ice(c,cptr%cps%snl(c)+1))
            else
              if ( cptr%cws%h2osno(c) > 0._rk8 ) then
                cptr%cps%snw_rds(c,0) = snw_rds_min
                cptr%cps%snw_rds(c,-nlevsno+1:-1) = 0._rk8
                cptr%cps%snw_rds_top(c) = spval
                cptr%cps%sno_liq_top(c) = spval
              else
                cptr%cps%snw_rds(c,:) = 0._rk8
                cptr%cps%snw_rds_top(c) = spval
                cptr%cps%sno_liq_top(c) = spval
              end if
            end if
          end do
        end if
      else
        call clm_readvar(ncid,'snw_rds',cptr%cps%snw_rds,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'snw_rds',cptr%cps%snw_rds,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - mss_bcpho

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'mss_bcpho', ['column','levsno'], &
              long_name='snow layer hydrophobic black carbon mass', &
              units='kg m-2', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'mss_bcpho') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          ! initial run, not restart: initialize mss_bcpho to zero
          do c = begc , endc
            cptr%cps%mss_bcpho(c,-nlevsno+1:0) = 0._rk8
          end do
        end if
      else
        call clm_readvar(ncid,'mss_bcpho',cptr%cps%mss_bcpho,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'mss_bcpho',cptr%cps%mss_bcpho,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - mss_bcphi

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'mss_bcphi', ['column','levsno'], &
              long_name='snow layer hydrophilic black carbon mass', &
              units='kg m-2', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'mss_bcphi') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          ! initial run, not restart: initialize mss_bcphi to zero
          do c = begc , endc
            cptr%cps%mss_bcphi(c,-nlevsno+1:0) = 0._rk8
          end do
        end if
      else
        call clm_readvar(ncid,'mss_bcphi',cptr%cps%mss_bcphi,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'mss_bcphi',cptr%cps%mss_bcphi,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - mss_ocpho

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'mss_ocpho', ['column','levsno'], &
              long_name='snow layer hydrophobic organic carbon mass', &
              units='kg m-2', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'mss_ocpho') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          ! initial run, not restart: initialize mss_ocpho to zero
          do c = begc , endc
            cptr%cps%mss_ocpho(c,-nlevsno+1:0) = 0._rk8
          end do
        end if
      else
        call clm_readvar(ncid,'mss_ocpho',cptr%cps%mss_ocpho,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'mss_ocpho',cptr%cps%mss_ocpho,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - mss_ocphi

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'mss_ocphi', ['column','levsno'], &
              long_name='snow layer hydrophilic organic carbon mass', &
              units='kg m-2', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'mss_ocphi') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          ! initial run, not restart: initialize mss_ocphi to zero
          do c = begc , endc
            cptr%cps%mss_ocphi(c,-nlevsno+1:0) = 0._rk8
          end do
        end if
      else
        call clm_readvar(ncid,'mss_ocphi',cptr%cps%mss_ocphi,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'mss_ocphi',cptr%cps%mss_ocphi,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - mss_dst1

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'mss_dst1', ['column','levsno'], &
              long_name='snow dust species 1 mass',units='kg m-2', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'mss_dst1') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          ! initial run, not restart: initialize mss_dst1 to zero
          do c = begc , endc
            cptr%cps%mss_dst1(c,-nlevsno+1:0) = 0._rk8
          end do
        end if
      else
        call clm_readvar(ncid,'mss_dst1',cptr%cps%mss_dst1,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'mss_dst1',cptr%cps%mss_dst1,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - mss_dst2

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'mss_dst2', ['column','levsno'], &
              long_name='snow dust species 2 mass',units='kg m-2', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'mss_dst2') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          ! initial run, not restart: initialize mss_dst2 to zero
          do c = begc , endc
            cptr%cps%mss_dst2(c,-nlevsno+1:0) = 0._rk8
          end do
        end if
      else
        call clm_readvar(ncid,'mss_dst2',cptr%cps%mss_dst2,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'mss_dst2',cptr%cps%mss_dst2,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - mss_dst3

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'mss_dst3', ['column','levsno'], &
              long_name='snow dust species 3 mass',units='kg m-2', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'mss_dst3') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          ! initial run, not restart: initialize mss_dst3 to zero
          do c = begc , endc
            cptr%cps%mss_dst3(c,-nlevsno+1:0) = 0._rk8
          end do
        end if
      else
        call clm_readvar(ncid,'mss_dst3',cptr%cps%mss_dst3,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'mss_dst3',cptr%cps%mss_dst3,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - mss_dst4

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'mss_dst4', ['column','levsno'], &
              long_name='snow dust species 4 mass',units='kg m-2', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'mss_dst4') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          ! initial run, not restart: initialize mss_dst4 to zero
          do c = begc , endc
            cptr%cps%mss_dst4(c,-nlevsno+1:0) = 0._rk8
          end do
        end if
      else
        call clm_readvar(ncid,'mss_dst4',cptr%cps%mss_dst4,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'mss_dst4',cptr%cps%mss_dst4,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - flx_absdv

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'flx_absdv', ['column ','levsno1'], &
              long_name='snow layer flux absorption factors (direct, VIS)', &
              units='1', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'flx_absdv') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          do_initsurfalb = .true.
        end if
      else
        call clm_readvar(ncid,'flx_absdv',cptr%cps%flx_absdv,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'flx_absdv',cptr%cps%flx_absdv,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - flx_absdn

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'flx_absdn', ['column ','levsno1'], &
              long_name='snow layer flux absorption factors (direct, NIR)', &
              units='1', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'flx_absdn') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          do_initsurfalb = .true.
        end if
      else
        call clm_readvar(ncid,'flx_absdn',cptr%cps%flx_absdn,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'flx_absdn',cptr%cps%flx_absdn,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - flx_absiv

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'flx_absiv', ['column ','levsno1'], &
              long_name='snow layer flux absorption factors (diffuse, VIS)', &
              units='1', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'flx_absiv') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          do_initsurfalb = .true.
        end if
      else
        call clm_readvar(ncid,'flx_absiv',cptr%cps%flx_absiv,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'flx_absiv',cptr%cps%flx_absiv,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - flx_absin

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'flx_absin', ['column ','levsno1'], &
              long_name='snow layer flux absorption factors (diffuse, NIR)', &
              units='1', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'flx_absin') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          do_initsurfalb = .true.
        end if
      else
        call clm_readvar(ncid,'flx_absin',cptr%cps%flx_absin,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'flx_absin',cptr%cps%flx_absin,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - albsnd_hst

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albsnd_hst',['column','numrad'], &
              long_name='snow albedo (direct)',units='1', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'albsnd_hst') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'albsnd_hst',cptr%cps%albsnd_hst,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albsnd_hst',cptr%cps%albsnd_hst,gcomm_column, switchdim=.true.)
    end if

    ! column type physical state variable - albsni_hst

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'albsni_hst',['column','numrad'], &
              long_name='snow albedo (diffuse)',units='1', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( lstart .and. .not. clm_check_var(ncid,'albsni_hst') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'albsni_hst',cptr%cps%albsni_hst,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'albsni_hst',cptr%cps%albsni_hst,gcomm_column, switchdim=.true.)
    end if

    ! column type water flux variable - qflx_snofrz_lyr

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'qflx_snofrz_lyr', &
              ['column','levsno'], &
              long_name='snow layer ice freezing rate',units='kg m-2 s-1', switchdim=.true.)
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'qflx_snofrz_lyr') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          do c = begc , endc
            cptr%cwf%qflx_snofrz_lyr(c,-nlevsno+1:0) = 0._rk8
          end do
        end if
      else
        call clm_readvar(ncid,'qflx_snofrz_lyr', &
                cptr%cwf%qflx_snofrz_lyr,gcomm_column, switchdim=.true.)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'qflx_snofrz_lyr', &
              cptr%cwf%qflx_snofrz_lyr,gcomm_column, switchdim=.true.)
    end if

    ! column type water flux variable - qflx_snow_melt

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'qflx_snow_melt',['column'], &
              long_name='net snow melt',units='mm s-1')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'qflx_snow_melt') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          ! initial run, not restart: initialize qflx_snow_melt to zero
          cptr%cwf%qflx_snow_melt = 0._rk8
        end if
      else
        call clm_readvar(ncid,'qflx_snow_melt', &
                cptr%cwf%qflx_snow_melt,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'qflx_snow_melt', &
              cptr%cwf%qflx_snow_melt,gcomm_column)
    end if

    ! gridcell type water flux variable - qflx_floodg

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'qflx_floodg',['gridcell'], &
              long_name='flood water flux',units='mm s-1')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'qflx_floodg') ) then
        if ( lstart ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          ! initial run, not restart: initialize flood to zero
          clm_a2l%forc_flood = 0._rk8
        end if
      else
        call clm_readvar(ncid,'qflx_floodg',clm_a2l%forc_flood,gcomm_gridcell)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'qflx_floodg',clm_a2l%forc_flood,gcomm_gridcell)
    end if

    ! initialize other variables that are derived from those
    ! stored in the restart buffer. (there may be a more appropriate
    ! place to do this, but functionally this works)

    if ( flag == 'read' ) then
      do j = -nlevsno+1 , 0
        do c = begc , endc
          ! mass concentrations of aerosols in snow
          if ( cptr%cws%h2osoi_ice(c,j)+cptr%cws%h2osoi_liq(c,j) > 0._rk8 ) then
            cptr%cps%mss_cnc_bcpho(c,j) = cptr%cps%mss_bcpho(c,j) / &
                    (cptr%cws%h2osoi_ice(c,j)+cptr%cws%h2osoi_liq(c,j))
            cptr%cps%mss_cnc_bcphi(c,j) = cptr%cps%mss_bcphi(c,j) / &
                    (cptr%cws%h2osoi_ice(c,j)+cptr%cws%h2osoi_liq(c,j))
            cptr%cps%mss_cnc_ocpho(c,j) = cptr%cps%mss_ocpho(c,j) / &
                    (cptr%cws%h2osoi_ice(c,j)+cptr%cws%h2osoi_liq(c,j))
            cptr%cps%mss_cnc_ocphi(c,j) = cptr%cps%mss_ocphi(c,j) / &
                    (cptr%cws%h2osoi_ice(c,j)+cptr%cws%h2osoi_liq(c,j))
            cptr%cps%mss_cnc_dst1(c,j) = cptr%cps%mss_dst1(c,j) / &
                    (cptr%cws%h2osoi_ice(c,j)+cptr%cws%h2osoi_liq(c,j))
            cptr%cps%mss_cnc_dst2(c,j) = cptr%cps%mss_dst2(c,j) / &
                    (cptr%cws%h2osoi_ice(c,j)+cptr%cws%h2osoi_liq(c,j))
            cptr%cps%mss_cnc_dst3(c,j) = cptr%cps%mss_dst3(c,j) / &
                    (cptr%cws%h2osoi_ice(c,j)+cptr%cws%h2osoi_liq(c,j))
            cptr%cps%mss_cnc_dst4(c,j) = cptr%cps%mss_dst4(c,j) / &
                    (cptr%cws%h2osoi_ice(c,j)+cptr%cws%h2osoi_liq(c,j))
          else
            cptr%cps%mss_cnc_bcpho(c,j) = 0._rk8
            cptr%cps%mss_cnc_bcphi(c,j) = 0._rk8
            cptr%cps%mss_cnc_ocpho(c,j) = 0._rk8
            cptr%cps%mss_cnc_ocphi(c,j) = 0._rk8
            cptr%cps%mss_cnc_dst1(c,j) = 0._rk8
            cptr%cps%mss_cnc_dst2(c,j) = 0._rk8
            cptr%cps%mss_cnc_dst3(c,j) = 0._rk8
            cptr%cps%mss_cnc_dst4(c,j) = 0._rk8
          end if
        end do
      end do
    end if
  end subroutine BiogeophysRest
  !
  ! Determine if the weights read in are exactly the same as those
  ! from surface dataset
  !
  logical function weights_exactly_the_same(pptr,wtgcell,wtlunit,wtcol)
    implicit none
    type(pft_type) , pointer :: pptr     ! pointer to pft derived subtype
    real(rk8) , intent(in) :: wtgcell(:) ! grid cell weights for each PFT
    real(rk8) , intent(in) :: wtlunit(:) ! land-unit weights for each PFT
    real(rk8) , intent(in) :: wtcol(:)   ! column weights for each PFT
    ! Check that weights are identical for all PFT's and all weight types
    if ( all( pptr%wtgcell(:) == wtgcell ) .and. &
         all( pptr%wtlunit(:) == wtlunit ) .and. &
         all( pptr%wtcol(:) == wtcol ) ) then
      weights_exactly_the_same = .true.
    else
      weights_exactly_the_same = .false.
    end if
  end function weights_exactly_the_same
  !
  ! Determine if the weights are within roundoff different from each other
  !
  logical function weights_within_roundoff_different(pptr,wtgcell,wtlunit,wtcol)
    implicit none
    type(pft_type) , pointer :: pptr ! pointer to pft derived subtype
    real(rk8) , intent(in) :: wtgcell(:) ! grid cell weights for each PFT
    real(rk8) , intent(in) :: wtlunit(:) ! land-unit weights for each PFT
    real(rk8) , intent(in) :: wtcol(:)   ! column weights for each PFT
    real(rk8) , parameter :: rndVal = 1.e-13_rk8
    ! If differences between all weights for each PFT and each weight type is
    ! less than or equal to double precision roundoff level
    ! weights are close
    if ( all(abs(pptr%wtgcell(:) - wtgcell) <= rndVal) .and. &
         all(abs(pptr%wtlunit(:) - wtlunit) <= rndVal) .and. &
         all(abs(pptr%wtcol(:)   - wtcol  ) <= rndVal) ) then
      weights_within_roundoff_different = .true.
    else
      weights_within_roundoff_different = .false.
    end if
  end function weights_within_roundoff_different
  !
  ! Determine if the weights read in are too different and should flag an error
  !
  logical function weights_tooDifferent(begp,endp,pptr,wtgcell,adiff,maxdiff)
    implicit none
    ! per-proc beginning and ending pft indices
    integer(ik4) , intent(in) :: begp, endp
    type(pft_type) , pointer :: pptr ! pointer to pft derived subtype
    ! grid cell weights for each PFT
    real(rk8) , intent(in) :: wtgcell(begp:endp)
    real(rk8) , intent(in) :: adiff ! tolerance of acceptible difference
    real(rk8) , intent(out) :: maxdiff ! maximum difference found
    integer(ik4) :: p     ! PFT index
    real(rk8) :: diff     ! difference in weights
    ! Assume weights are NOT different and only change if find
    ! weights too different
    weights_tooDifferent = .false.
    maxdiff = 0.0_rk8
    do p = begp , endp
      diff = abs(pptr%wtgcell(p) - wtgcell(p))
      if ( diff > maxdiff ) maxdiff = diff
      if ( diff > adiff ) weights_tooDifferent = .true.
    end do
  end function weights_tooDifferent

end module mod_clm_biogeophysrest
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
