module mod_clm_urbaninput

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: UrbanInputMod
! 
! !DESCRIPTION: 
! Read in input urban data - fill in data structure urbinp
!
! !USES:
  use mod_realkinds
  use mod_mpmessage
  use mod_nchelper
  use mod_dynparam
  use mod_mppparam
!
! !PUBLIC TYPES:
  implicit none
  save

  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: UrbanInput         ! Read in urban input data

  type urbinp_t
     real(rk8), pointer :: canyon_hwr(:,:)  
     real(rk8), pointer :: wtlunit_roof(:,:)  
     real(rk8), pointer :: wtroad_perv(:,:)  
     real(rk8), pointer :: em_roof(:,:)   
     real(rk8), pointer :: em_improad(:,:)  
     real(rk8), pointer :: em_perroad(:,:)  
     real(rk8), pointer :: em_wall(:,:)  
     real(rk8), pointer :: alb_roof_dir(:,:,:)  
     real(rk8), pointer :: alb_roof_dif(:,:,:)  
     real(rk8), pointer :: alb_improad_dir(:,:,:)  
     real(rk8), pointer :: alb_improad_dif(:,:,:)  
     real(rk8), pointer :: alb_perroad_dir(:,:,:)  
     real(rk8), pointer :: alb_perroad_dif(:,:,:)  
     real(rk8), pointer :: alb_wall_dir(:,:,:)  
     real(rk8), pointer :: alb_wall_dif(:,:,:)  
     real(rk8), pointer :: ht_roof(:,:)
     real(rk8), pointer :: wind_hgt_canyon(:,:)
     real(rk8), pointer :: tk_wall(:,:,:)
     real(rk8), pointer :: tk_roof(:,:,:)
     real(rk8), pointer :: tk_improad(:,:,:)
     real(rk8), pointer :: cv_wall(:,:,:)
     real(rk8), pointer :: cv_roof(:,:,:)
     real(rk8), pointer :: cv_improad(:,:,:)
     real(rk8), pointer :: thick_wall(:,:)
     real(rk8), pointer :: thick_roof(:,:)
     integer,  pointer :: nlev_improad(:,:)
     real(rk8), pointer :: t_building_min(:,:)
     real(rk8), pointer :: t_building_max(:,:)
  end type urbinp_t
  public urbinp_t

  type (urbinp_t)   , public :: urbinp        ! urban input derived type
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: UrbanInput
!
! !INTERFACE:
  subroutine UrbanInput(mode)
!
! !DESCRIPTION: 
! Allocate memory and read in urban input data
!
! !USES:
    use mod_clm_varpar, only : numrad, nlevurb, numurbl
    use mod_clm_varctl, only : fsurdat, single_column
    use mod_clm_type   , only : grlnd
    use mod_clm_decomp , only : get_proc_bounds
    use mod_clm_domain , only : ldomain
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: mode
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein July 2004
! Revised by Keith Oleson for netcdf input Jan 2008
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: ncid                  ! netcdf id
    integer :: begg,endg             ! start/stop gridcells
    integer :: nw,n,k,i,j,ni,nj,ns   ! indices
    integer :: nlevurb_i             ! input grid: number of urban vertical levels
    integer :: numrad_i              ! input grid: number of solar bands (VIS/NIR)
    integer :: numurbl_i             ! input grid: number of urban landunits
    integer :: ier,ret               ! error status
    logical :: isgrid2d              ! true => file is 2d 
    logical :: readvar               ! true => variable is on dataset
    logical :: has_numurbl           ! true => numurbl dimension is on dataset
    logical :: has_nsolar            ! true => nsolar dimension is on dataset
    character(len=32) :: subname = 'UrbanInput' ! subroutine name
!-----------------------------------------------------------------------

    if ( nlevurb == 0 ) return

    call get_proc_bounds(begg,endg)

    if (mode == 'initialize') then

       ! Read urban data
       
       if (myid == italk) then
          write(stdout,*)' Reading in urban input data from fsurdat file ...'
       end if
       
       call openfile_withname(fsurdat,ncid)

       if (myid == italk) then
          write(stdout,*) subname,trim(fsurdat)
       end if

       call ncd_inqdim(ncid,'nsolar',lexist=has_nsolar)

       ! Check whether this file has new-format urban data
       call ncd_inqdim(ncid,'numurbl',lexist=has_numurbl)

       ! If file doesn't have numurbl, then it is old-format urban;
       ! in this case, set nlevurb to zero
       if (.not. has_numurbl) then
         nlevurb = 0
         write(stdout,*)'PCT_URBAN is not multi-density, nlevurb set to 0'
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
       call check_dims(ncid,ni,nj)
       ns = ni*nj
       if (ldomain%ns /= ns .or. ldomain%ni /= ni .or. ldomain%nj /= nj) then
          write(stderr,*)trim(subname), &
                    'ldomain and input file do not match dims '
          write(stderr,*)trim(subname), 'ldomain%ni,ni,= ',ldomain%ni,ni
          write(stderr,*)trim(subname), 'ldomain%nj,nj,= ',ldomain%nj,nj
          write(stderr,*)trim(subname), 'ldomain%ns,ns,= ',ldomain%ns,ns
          call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
       call ncd_inqdim(ncid,'nlevurb',dlen=nlevurb_i)
       if (nlevurb_i /= nlevurb) then
          write(stderr,*)trim(subname)// ': parameter nlevurb= ',nlevurb, &
               'does not equal input dataset nlevurb= ',nlevurb_i
          call fatal(__FILE__,__LINE__,'clm now stopping')
       endif
       call ncd_inqdim(ncid,'numrad',dlen=numrad_i)
       if (numrad_i /= numrad) then
          write(stderr,*)trim(subname)// ': parameter numrad= ',numrad, &
               'does not equal input dataset numrad= ',numrad_i
          call fatal(__FILE__,__LINE__,'clm now stopping')
       endif
       call ncd_inqdim(ncid,'numurbl',dlen=numurbl_i)
       if (numurbl_i /= numurbl) then
          write(stderr,*)trim(subname)// ': parameter numurbl= ',numurbl, &
               'does not equal input dataset numurbl= ',numurbl_i
          call fatal(__FILE__,__LINE__,'clm now stopping')
       endif

       call read_var2d_static(ncid,'CANYON_HWR',ni,nj,urbinp%canyon_hwr)
       call read_var2d_static(ncid,'WTLUNIT_ROOF',ni,nj,urbinp%wtlunit_roof)
       call read_var2d_static(ncid,'WTROAD_PERV',ni,nj,urbinp%wtroad_perv)
       call read_var2d_static(ncid,'EM_ROOF',ni,nj,urbinp%em_roof)
       call read_var2d_static(ncid,'EM_IMPROAD',ni,nj,urbinp%em_improad)
       call read_var2d_static(ncid,'EM_PERROAD',ni,nj,urbinp%em_perroad)
       call read_var2d_static(ncid,'EM_WALL',ni,nj,urbinp%em_wall)
       call read_var2d_static(ncid,'HT_ROOF',ni,nj,urbinp%ht_roof)
       call read_var2d_static(ncid,'WIND_HGT_CANYON',ni,nj,urbinp%wind_hgt_canyon)
       call read_var2d_static(ncid,'THICK_WALL',ni,nj,urbinp%thick_wall)
       call read_var2d_static(ncid,'THICK_ROOF',ni,nj,urbinp%thick_roof)
       call read_var2d_static(ncid,'NLEV_IMPROAD',ni,nj,urbinp%nlev_improad)
       call read_var2d_static(ncid,'T_BUILDING_MIN',ni,nj,urbinp%t_building_min)
       call read_var2d_static(ncid,'T_BUILDING_MAX',ni,nj,urbinp%t_building_max)

       if ( has_nsolar ) then
         call fatal(__FILE__,__LINE__, &
             trim(subname)//' ERROR: Cannot handle nsolar here' )
         ! not implemented
       else
         call read_var3d_static(ncid,'ALB_IMPROAD_DIR',ni,nj,numrad,urbinp%alb_improad_dir)
         call read_var3d_static(ncid,'ALB_IMPROAD_DIF',ni,nj,numrad,urbinp%alb_improad_dif)
         call read_var3d_static(ncid,'ALB_PERROAD_DIR',ni,nj,numrad,urbinp%alb_perroad_dir)
         call read_var3d_static(ncid,'ALB_PERROAD_DIF',ni,nj,numrad,urbinp%alb_perroad_dif)
         call read_var3d_static(ncid,'ALB_ROOF_DIR',ni,nj,numrad,urbinp%alb_roof_dir)
         call read_var3d_static(ncid,'ALB_ROOF_DIF',ni,nj,numrad,urbinp%alb_roof_dif)
         call read_var3d_static(ncid,'ALB_WALL_DIR',ni,nj,numrad,urbinp%alb_wall_dir)
         call read_var3d_static(ncid,'ALB_WALL_DIF',ni,nj,numrad,urbinp%alb_wall_dif)
       end if
       call read_var3d_static(ncid,'TK_IMPROAD',ni,nj,nlevurb,urbinp%tk_improad)
       call read_var3d_static(ncid,'TK_ROOF',ni,nj,nlevurb,urbinp%tk_roof)
       call read_var3d_static(ncid,'TK_WALL',ni,nj,nlevurb,urbinp%tk_wall)
       call read_var3d_static(ncid,'CV_IMPROAD',ni,nj,nlevurb,urbinp%cv_improad)
       call read_var3d_static(ncid,'CV_ROOF',ni,nj,nlevurb,urbinp%cv_roof)
       call read_var3d_static(ncid,'CV_WALL',ni,nj,nlevurb,urbinp%cv_wall)

       call closefile(ncid)
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
