module mod_clm_urbaninput
  !
  ! Read in input urban data - fill in data structure urbinp
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_dynparam
  use mod_mppparam
  use mod_clm_nchelper
  use mod_clm_decomp , only : get_proc_bounds , gcomm_gridcell

  implicit none

  private

  save
  public :: UrbanInput ! Read in urban input data

  type urbinp_t
    real(rk8) , pointer , dimension(:,:) :: canyon_hwr
    real(rk8) , pointer , dimension(:,:) :: wtlunit_roof
    real(rk8) , pointer , dimension(:,:) :: wtroad_perv
    real(rk8) , pointer , dimension(:,:) :: em_roof
    real(rk8) , pointer , dimension(:,:) :: em_improad
    real(rk8) , pointer , dimension(:,:) :: em_perroad
    real(rk8) , pointer , dimension(:,:) :: em_wall
    real(rk8) , pointer , dimension(:,:,:) :: alb_roof_dir
    real(rk8) , pointer , dimension(:,:,:) :: alb_roof_dif
    real(rk8) , pointer , dimension(:,:,:) :: alb_improad_dir
    real(rk8) , pointer , dimension(:,:,:) :: alb_improad_dif
    real(rk8) , pointer , dimension(:,:,:) :: alb_perroad_dir
    real(rk8) , pointer , dimension(:,:,:) :: alb_perroad_dif
    real(rk8) , pointer , dimension(:,:,:) :: alb_wall_dir
    real(rk8) , pointer , dimension(:,:,:) :: alb_wall_dif
    real(rk8) , pointer , dimension(:,:) :: ht_roof
    real(rk8) , pointer , dimension(:,:) :: wind_hgt_canyon
    real(rk8) , pointer , dimension(:,:,:) :: tk_wall
    real(rk8) , pointer , dimension(:,:,:) :: tk_roof
    real(rk8) , pointer , dimension(:,:,:) :: tk_improad
    real(rk8) , pointer , dimension(:,:,:) :: cv_wall
    real(rk8) , pointer , dimension(:,:,:) :: cv_roof
    real(rk8) , pointer , dimension(:,:,:) :: cv_improad
    real(rk8) , pointer , dimension(:,:) :: thick_wall
    real(rk8) , pointer , dimension(:,:) :: thick_roof
    integer(ik4) , pointer , dimension(:,:) :: nlev_improad
    real(rk8) , pointer , dimension(:,:) :: t_building_min
    real(rk8) , pointer , dimension(:,:) :: t_building_max
  end type urbinp_t

  public urbinp_t

  type(urbinp_t) , public :: urbinp ! urban input derived type

  contains
  !
  ! Allocate memory and read in urban input data
  !
  subroutine UrbanInput(mode)
    use mod_clm_varpar , only : numrad , nlevurb , numurbl
    use mod_clm_varctl , only : fsurdat
    use mod_clm_domain , only : ldomain
    implicit none
    character(len=*) , intent(in) :: mode
    type(clm_filetype) :: ncid    ! netcdf id
    integer(ik4) :: begg , endg   ! start/stop gridcells
    integer(ik4) :: nlevurb_i ! input grid: number of urban vertical levels
    integer(ik4) :: numrad_i  ! input grid: number of solar bands (VIS/NIR)
    integer(ik4) :: numurbl_i ! input grid: number of urban landunits
    integer(ik4) :: ier       ! error status
    logical :: has_numurbl ! true => numurbl dimension is on dataset
    logical :: has_nsolar  ! true => nsolar dimension is on dataset
    character(len=32) :: subname = 'UrbanInput' ! subroutine name

    if ( nlevurb == 0 ) return

    call get_proc_bounds(begg,endg)

    if (mode == 'initialize') then

      ! Read urban data

      if (myid == italk) then
        write(stdout,*) ' Reading in urban input data from surface data file ...'
      end if

      call clm_openfile(fsurdat,ncid)

      if (myid == italk) then
        write(stdout,*) ' ',trim(subname),' : ',trim(fsurdat)
      end if

      ! Check whether this file has new-format urban data
      has_nsolar = clm_check_dim(ncid,'nsolar')
      has_numurbl = clm_check_dim(ncid,'numurbl')

      ! If file doesn't have numurbl, then it is old-format urban;
      ! in this case, set nlevurb to zero
      if ( .not. enable_urban_landunit ) then
        nlevurb = 0
        if ( myid == italk ) then
          write(stdout,*) 'Disabled urban area.'
        end if
      end if
      if (.not. has_numurbl) then
        nlevurb = 0
        if ( myid == italk ) then
          write(stdout,*) 'PCT_URBAN is not multi-density, nlevurb set to 0'
        end if
      end if

      if ( nlevurb == 0 ) return

      ! Allocate dynamic memory
      allocate(urbinp%canyon_hwr(begg:endg,numurbl), &
               urbinp%wtlunit_roof(begg:endg,numurbl), &
               urbinp%wtroad_perv(begg:endg,numurbl), &
               urbinp%em_roof(begg:endg,numurbl), &
               urbinp%em_improad(begg:endg,numurbl), &
               urbinp%em_perroad(begg:endg,numurbl), &
               urbinp%em_wall(begg:endg,numurbl), &
               urbinp%alb_roof_dir(begg:endg,numurbl,numrad), &
               urbinp%alb_roof_dif(begg:endg,numurbl,numrad), &
               urbinp%alb_improad_dir(begg:endg,numurbl,numrad), &
               urbinp%alb_perroad_dir(begg:endg,numurbl,numrad), &
               urbinp%alb_improad_dif(begg:endg,numurbl,numrad), &
               urbinp%alb_perroad_dif(begg:endg,numurbl,numrad), &
               urbinp%alb_wall_dir(begg:endg,numurbl,numrad), &
               urbinp%alb_wall_dif(begg:endg,numurbl,numrad), &
               urbinp%ht_roof(begg:endg,numurbl), &
               urbinp%wind_hgt_canyon(begg:endg,numurbl), &
               urbinp%tk_wall(begg:endg,numurbl,nlevurb), &
               urbinp%tk_roof(begg:endg,numurbl,nlevurb), &
               urbinp%tk_improad(begg:endg,numurbl,nlevurb), &
               urbinp%cv_wall(begg:endg,numurbl,nlevurb), &
               urbinp%cv_roof(begg:endg,numurbl,nlevurb), &
               urbinp%cv_improad(begg:endg,numurbl,nlevurb), &
               urbinp%thick_wall(begg:endg,numurbl), &
               urbinp%thick_roof(begg:endg,numurbl), &
               urbinp%nlev_improad(begg:endg,numurbl), &
               urbinp%t_building_min(begg:endg,numurbl), &
               urbinp%t_building_max(begg:endg,numurbl), &
               stat=ier)
      if (ier /= 0) then
        write(stderr,*)'initUrbanInput: allocation error '
        call fatal(__FILE__,__LINE__,'clm now stopping')
      endif
      call clm_inqdim(ncid,'nlevurb',nlevurb_i)
      if (nlevurb_i /= nlevurb) then
        write(stderr,*)trim(subname)// ': parameter nlevurb= ',nlevurb, &
             'does not equal input dataset nlevurb= ',nlevurb_i
        call fatal(__FILE__,__LINE__,'clm now stopping')
      endif
      call clm_inqdim(ncid,'numrad',numrad_i)
      if (numrad_i /= numrad) then
        write(stderr,*)trim(subname)// ': parameter numrad= ',numrad, &
             'does not equal input dataset numrad= ',numrad_i
        call fatal(__FILE__,__LINE__,'clm now stopping')
      endif
      call clm_inqdim(ncid,'numurbl',numurbl_i)
      if (numurbl_i /= numurbl) then
        write(stderr,*)trim(subname)// ': parameter numurbl= ',numurbl, &
             'does not equal input dataset numurbl= ',numurbl_i
        call fatal(__FILE__,__LINE__,'clm now stopping')
      endif

      call clm_readvar(ncid,'CANYON_HWR',urbinp%canyon_hwr,gcomm_gridcell)
      call clm_readvar(ncid,'WTLUNIT_ROOF',urbinp%wtlunit_roof,gcomm_gridcell)
      call clm_readvar(ncid,'WTROAD_PERV',urbinp%wtroad_perv,gcomm_gridcell)
      call clm_readvar(ncid,'EM_ROOF',urbinp%em_roof,gcomm_gridcell)
      call clm_readvar(ncid,'EM_IMPROAD',urbinp%em_improad,gcomm_gridcell)
      call clm_readvar(ncid,'EM_PERROAD',urbinp%em_perroad,gcomm_gridcell)
      call clm_readvar(ncid,'EM_WALL',urbinp%em_wall,gcomm_gridcell)
      call clm_readvar(ncid,'HT_ROOF',urbinp%ht_roof,gcomm_gridcell)
      call clm_readvar(ncid,'WIND_HGT_CANYON', &
        urbinp%wind_hgt_canyon,gcomm_gridcell)
      call clm_readvar(ncid,'THICK_WALL',urbinp%thick_wall,gcomm_gridcell)
      call clm_readvar(ncid,'THICK_ROOF',urbinp%thick_roof,gcomm_gridcell)
      call clm_readvar(ncid,'NLEV_IMPROAD',urbinp%nlev_improad,gcomm_gridcell)
      call clm_readvar(ncid,'T_BUILDING_MIN', &
        urbinp%t_building_min,gcomm_gridcell)
      call clm_readvar(ncid,'T_BUILDING_MAX', &
        urbinp%t_building_max,gcomm_gridcell)

      if ( has_nsolar ) then
        call fatal(__FILE__,__LINE__, &
            trim(subname)//' ERROR: Cannot handle nsolar here' )
        ! not implemented
      else
        call clm_readvar(ncid,'ALB_IMPROAD_DIR', &
          urbinp%alb_improad_dir,gcomm_gridcell)
        call clm_readvar(ncid,'ALB_IMPROAD_DIF', &
          urbinp%alb_improad_dif,gcomm_gridcell)
        call clm_readvar(ncid,'ALB_PERROAD_DIR', &
          urbinp%alb_perroad_dir,gcomm_gridcell)
        call clm_readvar(ncid,'ALB_PERROAD_DIF', &
          urbinp%alb_perroad_dif,gcomm_gridcell)
        call clm_readvar(ncid,'ALB_ROOF_DIR', &
          urbinp%alb_roof_dir,gcomm_gridcell)
        call clm_readvar(ncid,'ALB_ROOF_DIF', &
          urbinp%alb_roof_dif,gcomm_gridcell)
        call clm_readvar(ncid,'ALB_WALL_DIR', &
          urbinp%alb_wall_dir,gcomm_gridcell)
        call clm_readvar(ncid,'ALB_WALL_DIF', &
          urbinp%alb_wall_dif,gcomm_gridcell)
      end if
      call clm_readvar(ncid,'TK_IMPROAD',urbinp%tk_improad,gcomm_gridcell)
      call clm_readvar(ncid,'TK_ROOF',urbinp%tk_roof,gcomm_gridcell)
      call clm_readvar(ncid,'TK_WALL',urbinp%tk_wall,gcomm_gridcell)
      call clm_readvar(ncid,'CV_IMPROAD',urbinp%cv_improad,gcomm_gridcell)
      call clm_readvar(ncid,'CV_ROOF',urbinp%cv_roof,gcomm_gridcell)
      call clm_readvar(ncid,'CV_WALL',urbinp%cv_wall,gcomm_gridcell)

      call clm_closefile(ncid)
      if (myid == italk) then
        write(stdout,*)' Sucessfully read urban input data'
        write(stdout,*)
      end if

    else if (mode == 'finalize') then

      if ( nlevurb == 0 ) return

      deallocate(urbinp%canyon_hwr, &
                 urbinp%wtlunit_roof, &
                 urbinp%wtroad_perv, &
                 urbinp%em_roof, &
                 urbinp%em_improad, &
                 urbinp%em_perroad, &
                 urbinp%em_wall, &
                 urbinp%alb_roof_dir, &
                 urbinp%alb_roof_dif, &
                 urbinp%alb_improad_dir, &
                 urbinp%alb_perroad_dir, &
                 urbinp%alb_improad_dif, &
                 urbinp%alb_perroad_dif, &
                 urbinp%alb_wall_dir, &
                 urbinp%alb_wall_dif, &
                 urbinp%ht_roof, &
                 urbinp%wind_hgt_canyon, &
                 urbinp%tk_wall, &
                 urbinp%tk_roof, &
                 urbinp%tk_improad, &
                 urbinp%cv_wall, &
                 urbinp%cv_roof, &
                 urbinp%cv_improad, &
                 urbinp%thick_wall, &
                 urbinp%thick_roof, &
                 urbinp%nlev_improad, &
                 urbinp%t_building_min, &
                 urbinp%t_building_max, &
                 stat=ier)
      if (ier /= 0) then
        write(stderr,*)'initUrbanInput: deallocation error '
        call fatal(__FILE__,__LINE__,'clm now stopping')
      endif

    else
      write(stderr,*)'initUrbanInput error: mode ',trim(mode),' not supported '
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

  end subroutine UrbanInput

end module mod_clm_urbaninput
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
