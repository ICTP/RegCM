module mod_clm_surfrd
  !
  ! Contains methods for reading in surface data file and determining
  ! two-dimensional subgrid weights as well as writing out new surface
  ! dataset. When reading in the surface dataset, determines array
  ! which sets the PFT for each of the [maxpatch] patches and
  ! array which sets the relative abundance of the PFT.
  ! Also fills in the PFTs for vegetated portion of each grid cell.
  ! Fractional areas for these points pertain to "vegetated"
  ! area not to total grid area. Need to adjust them for fraction of grid
  ! that is vegetated. Also fills in urban, lake, wetland, and glacier patches.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_memutil
  use mod_mpmessage
  use mod_dynparam , only : myid
#if (defined CNDV)
  use mod_dynparam , only : enable_dv_baresoil
#endif
  use mod_mppparam
  use mod_clm_nchelper
  use mod_clm_varpar , only : nlevsoifl , numpft , maxpatch_pft , &
         maxpatch , npatch_urban_tbd , npatch_urban_hd , npatch_urban_md , &
         numurbl , npatch_lake , npatch_wet , npatch_glacier , maxpatch_urb
  use mod_clm_varsur , only : wtxy , vegxy , pctspec
  use mod_clm_decomp , only : get_proc_bounds , procinfo , numg
  use mod_clm_decomp , only : gcomm_gridcell
  use mod_clm_type

  implicit none

  private

  save

  ! Read grid/ladnfrac data into domain (after domain decomp)
  public :: surfrd_get_grid
  ! Read surface dataset and determine subgrid weights
  public :: surfrd_get_data

  ! If prognostic crops is turned on
  logical , public :: crop_prog = .false.

  private :: surfrd_wtxy_special
  private :: surfrd_wtxy_veg_all
  private :: surfrd_wtxy_veg_dgvm

  ! default multiplication factor for epsilon for error checks
  real(rk8) , private, parameter :: eps_fact = 2._rk8

  contains
  !
  ! THIS IS CALLED AFTER THE DOMAIN DECOMPOSITION HAS BEEN CREATED
  ! Read the surface dataset grid related information:
  ! o real latitude  of grid cell (degrees)
  ! o real longitude of grid cell (degrees)
  !
  subroutine surfrd_get_grid(ldomain, filename)
    use mod_clm_varcon , only : spval , re
    use mod_clm_domain , only : domain_type , domain_init , domain_clean
    use mod_clm_decomp , only : get_proc_bounds
    use mod_clm_atmlnd , only : adomain
    implicit none
    type(domain_type) , intent(inout) :: ldomain   ! domain to init
    character(len=*) , intent(in)    :: filename  ! grid filename
    type(clm_filetype) :: ncid     ! netcdf id
    integer(ik4) :: ibeg           ! local beg index
    integer(ik4) :: iend           ! local end index
    integer(ik4) :: ni , nj , ns   ! size of grid on file
    integer(ik4) :: inni , innj    ! size of grid on file
    real(rk4) , allocatable , dimension(:,:) :: mask
    character(len=32) :: subname = 'surfrd_get_grid'     ! subroutine name

    if (myid == italk) then
      if (filename == ' ') then
        write(stderr,*) trim(subname),' ERROR: filename must be specified '
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
    end if

    call clm_openfile(filename,ncid)

    ! Determine dimensions

    call clm_inqdim(ncid,'jx',inni)
    call clm_inqdim(ncid,'iy',innj)

    allocate(mask(inni,innj))
    call clm_readvar(ncid,'mask',mask)

    call clm_closefile(ncid)

    ni = inni - 3
    nj = innj - 3
    ns = ni*nj

    call get_proc_bounds(ibeg,iend)
    call domain_init(ldomain,ni=ni,nj=nj,nbeg=ibeg,nend=iend)

    allocate(procinfo%gcmask(ni,nj))
    procinfo%gcmask(:,:) = (mask(2:inni-2,2:innj-2) > 0.0)
    deallocate(mask)

    ! We receive only GOOD land points !

    ldomain%frac = 100.0_rk8
    ldomain%mask = 1
    ldomain%lonc = adomain%xlon
    ldomain%latc = adomain%xlat
    ldomain%topo = adomain%topo
    ldomain%area = adomain%area*.000001_rk8 ! square kilometers
    ldomain%pftm = ldomain%mask

    ! Check lat limited to -90,90

    if (minval(ldomain%latc) < -90.0_rk8 .or. &
        maxval(ldomain%latc) >  90.0_rk8) then
      write(stderr,*) trim(subname),' WARNING: lat/lon min/max is ', &
           minval(ldomain%latc),maxval(ldomain%latc)
      call fatal(__FILE__,__LINE__, &
          trim(subname)//' ERROR: lat is outside [-90,90]' )
    end if

  end subroutine surfrd_get_grid
  !
  ! Read the surface dataset and create subgrid weights.
  ! The model's surface dataset recognizes 6 basic land cover types within
  ! a grid cell: lake, wetland, urban, glacier, glacier_mec and vegetated.
  ! The vegetated portion of the grid cell is comprised of up to [maxpatch_pft]
  ! PFTs. These subgrid patches are read in explicitly for each grid cell.
  ! This is in contrast to LSMv1, where the PFTs were built implicitly from
  ! biome types.
  !    o real latitude  of grid cell (degrees)
  !    o real longitude of grid cell (degrees)
  !    o integer(ik4) surface type: 0 = ocean or 1 = land
  !    o integer(ik4) soil color (1 to 20) for use with soil albedos
  !    o real soil texture, %sand, for thermal and hydraulic properties
  !    o real soil texture, %clay, for thermal and hydraulic properties
  !    o real % of cell covered by lake    for use as subgrid patch
  !    o real % of cell covered by wetland for use as subgrid patch
  !    o real % of cell that is urban      for use as subgrid patch
  !    o real % of cell that is glacier    for use as subgrid patch
  !    o real % of cell that is glacier_mec for use as subgrid patch
  !    o integer(ik4) PFTs
  !    o real % abundance PFTs (as a percent of vegetated area)
  subroutine surfrd_get_data(ldomain, lfsurdat)
    use mod_clm_varctl , only : allocate_all_vegpfts
#if (defined CNDV)
    use mod_clm_varctl , only : create_crop_landunit
#endif
    use mod_clm_pftvarcon , only : noveg
    use mod_clm_domain , only : domain_type , domain_init , domain_clean
    implicit none
    ! land domain associated with wtxy
    type(domain_type) , intent(inout) :: ldomain
    character(len=*) , intent(inout) :: lfsurdat ! surface dataset filename
    ! local domain associated with surface dataset
    integer(ik4) :: n , ierr                 ! loop indices
    type(clm_filetype) :: ncid               ! netcdf id
    integer(ik4) :: begg , endg              ! beg,end gridcell indices
    real(rkx) , allocatable , dimension(:) :: xclon , yclat
    character(len=32) :: subname = 'surfrd_get_data'    ! subroutine name

    if (myid == italk) then
      write(stdout,*) 'Attempting to read surface boundary data .....'
      if (lfsurdat == ' ') then
        write(stderr,*)'lfsurdat must be specified'
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
    end if

    call get_proc_bounds(begg,endg)
    allocate(pctspec(begg:endg))

    allocate(yclat(begg:endg))
    allocate(xclon(begg:endg))

    vegxy(:,:) = noveg
    wtxy(:,:)  = 0._rk8
    pctspec(:) = 0._rk8

    ! Read surface data

    call clm_openfile(lfsurdat,ncid)

    call clm_readvar(ncid,'xclon',xclon,gcomm_gridcell)
    call clm_readvar(ncid,'yclat',yclat,gcomm_gridcell)

    ierr = 0
    do n = begg , endg
      if ( abs(real(xclon(n),rk4)-real(ldomain%lonc(n),rk4)) > 10.e-7_rk8 .or. &
           abs(real(yclat(n),rk4)-real(ldomain%latc(n),rk4)) > 10.e-7_rk8 ) then
        write(stderr,*) 'ERROR coordinates at n ', &
            n, real(xclon(n),rk4), real(yclat(n),rk4) , &
               real(ldomain%lonc(n),rk4), real(ldomain%latc(n),rk4)
        ierr = 1
      end if
    end do
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'clm now stopping')
    else
      if ( myid == italk ) then
        write(stdout,*) 'Checked LAT/LON compliance'
      end if
    end if

    ! Obtain special landunit info

    call surfrd_wtxy_special(ncid,ldomain)

    ! Obtain vegetated landunit info

#if (defined CNDV)
    if ( enable_dv_baresoil ) then
      ! Reset to bare ground the vegetation
      call surfrd_wtxy_veg_dgvm()
    else
      call surfrd_wtxy_veg_all(ncid,ldomain)
    end if
#else
    if (allocate_all_vegpfts) then
      call surfrd_wtxy_veg_all(ncid,ldomain)
    else
      call fatal(__FILE__,__LINE__, &
              trim(subname) // 'only allocate_all_vegpfts is supported')
    end if
#endif
    call clm_closefile(ncid)

    if ( myid == italk )then
      write(stdout,*) 'Successfully read surface boundary data'
      write(stdout,*)
    end if

  end subroutine surfrd_get_data
  !
  ! Determine weight with respect to gridcell of all special "pfts" as well
  ! as soil color and percent sand and clay
  !
  subroutine surfrd_wtxy_special(ncid,ldomain)
    use mod_clm_pftvarcon , only : noveg
    use mod_clm_urbaninput , only : urbinp
    use mod_clm_varpar , only : nlevurb
    use mod_clm_varcon , only : udens_base , udens_tbd , udens_hd , udens_md
    use mod_clm_domain , only : domain_type
    implicit none
    type(clm_filetype) , intent(inout) :: ncid  ! netcdf id
    type(domain_type) , intent(inout) :: ldomain
    integer(ik4)  :: n , nl , nurb   ! indices
    integer(ik4)  :: begg , endg     ! gcell beg/end
    logical :: found                 ! temporary for error check
    integer(ik4) :: nindx            ! temporary for error check
    integer(ik4) :: nlev             ! level
    integer(ik4) :: dindx            ! temporary for error check
    ! percent of grid cell is glacier
    real(rk8) , pointer , dimension(:) :: pctgla
    ! percent of grid cell is lake
    real(rk8) , pointer , dimension(:) :: pctlak
    ! percent of grid cell is wetland
    real(rk8) , pointer , dimension(:) :: pctwet
    ! percent of grid cell is urbanized
    real(rk8) , pointer , dimension(:,:) :: pcturb
    ! percent of grid cell is urban (sum over density classes)
    real(rk8) , pointer , dimension(:) :: pcturb_tot
    character(len=32) :: subname = 'surfrd_wtxy_special'  ! subroutine name

    call get_proc_bounds(begg,endg)

    allocate(pctgla(begg:endg),pctlak(begg:endg),pctwet(begg:endg))
    allocate(pcturb(begg:endg,numurbl),pcturb_tot(begg:endg))

    if ( .not. clm_check_dimlen(ncid,'nlevsoi',nlevsoifl) ) then
      call fatal(__FILE__,__LINE__, &
        trim(subname)//' ERROR: nlevsoi not matching nlevsoifl' )
    end if

    ! Obtain non-grid surface properties of surface dataset other
    ! than percent pft

    if ( clm_check_var(ncid,'PCT_WETLAND') ) then
      call clm_readvar(ncid,'PCT_WETLAND',pctwet,gcomm_gridcell)
    else
      call fatal(__FILE__,__LINE__, &
        trim(subname)//' ERROR: PCT_WETLAND  NOT on surfdata file' )
    end if

    if ( clm_check_var(ncid,'PCT_LAKE') ) then
      call clm_readvar(ncid,'PCT_LAKE',pctlak,gcomm_gridcell)
    else
      call fatal(__FILE__,__LINE__, &
        trim(subname)//' ERROR: PCT_LAKE NOT on surfdata file' )
    end if

    if ( clm_check_var(ncid,'PCT_GLACIER') ) then
      call clm_readvar(ncid,'PCT_GLACIER',pctgla,gcomm_gridcell)
    else
      call fatal(__FILE__,__LINE__, &
        trim(subname)//' ERROR: PCT_GLACIER NOT on surfdata file' )
    end if

    ! If PCT_URBAN is not multi-density then set pcturb and nlevurb to zero
    if (nlevurb == 0) then
      pcturb = 0._rk8
      if ( myid == italk ) then
        write(stdout,*)'PCT_URBAN is not multi-density, pcturb set to 0'
      end if
    else
      if ( clm_check_var(ncid,'PCT_URBAN') ) then
        call clm_readvar(ncid,'PCT_URBAN',pcturb,gcomm_gridcell)
      else
        call fatal(__FILE__,__LINE__, &
          trim(subname)//' ERROR: PCT_URBAN NOT on surfdata file' )
      end if
    end if
    if ( nlevurb == 0 )then
      if ( any(pcturb > 0.0_rk8) ) then
        call fatal(__FILE__,__LINE__, &
          trim(subname)//' ERROR: PCT_URBAN MUST be zero when nlevurb=0' )
      end if
    end if

    pcturb_tot(:) = 0._rk8
    do n = 1 , numurbl
      do nl = begg , endg
        pcturb_tot(nl) = pcturb_tot(nl) + pcturb(nl,n)
      end do
    end do

    where ( pctwet < 1.0e-4_rk8 )
      pctwet = 0.0_rk8
    end where
    where ( pctlak < 1.0e-4_rk8 )
      pctlak = 0.0_rk8
    end where
    where ( pctgla < 1.0e-4_rk8 )
      pctgla = 0.0_rk8
    end where
    where ( pcturb_tot < 1.0e-4_rk8 )
      pcturb_tot = 0.0_rk8
    end where

    pctspec = pctwet + pctlak + pcturb_tot + pctgla

    ! Error check: glacier, lake, wetland, urban sum must be less than 100

    found = .false.
    do nl = begg , endg
      if (pctspec(nl) > 100.0001_rk8 ) then
        found = .true.
        nindx = nl
        exit
      end if
      if (found) exit
    end do
    if ( found ) then
      write(stderr,*) 'surfrd error: SPECIAL classes > 100 for nl=',nindx
      write(stderr,*) 'pctwet = ', pctwet(nl)
      write(stderr,*) 'pctlak = ', pctlak(nl)
      write(stderr,*) 'pctgla = ', pctgla(nl)
      write(stderr,*) 'pcturb = ', pcturb_tot(nl)
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! Determine veg and wtxy for special landunits

    do nl = begg , endg
      vegxy(nl,npatch_lake)   = noveg
      wtxy(nl,npatch_lake)    = pctlak(nl)/100._rk8

      vegxy(nl,npatch_wet)    = noveg
      wtxy(nl,npatch_wet)     = pctwet(nl)/100._rk8

      vegxy(nl,npatch_glacier)= noveg
      wtxy(nl,npatch_glacier) = pctgla(nl)/100._rk8

      ! Initialize urban tall building district weights
      n = udens_tbd - udens_base
      do nurb = npatch_urban_tbd , npatch_urban_hd-1
        vegxy(nl,nurb) = noveg
        wtxy(nl,nurb)  = pcturb(nl,n) / 100._rk8
      end do
      if ( pcturb(nl,n) > 0.0_rk8 )then
        wtxy(nl,npatch_urban_tbd) = wtxy(nl,npatch_urban_tbd) * &
                urbinp%wtlunit_roof(nl,n)
        wtxy(nl,npatch_urban_tbd+1) = wtxy(nl,npatch_urban_tbd+1) * &
                (1.0_rk8 - urbinp%wtlunit_roof(nl,n))/3.0_rk8
        wtxy(nl,npatch_urban_tbd+2) = wtxy(nl,npatch_urban_tbd+2) * &
                (1.0_rk8 - urbinp%wtlunit_roof(nl,n))/3.0_rk8
        wtxy(nl,npatch_urban_tbd+3) = wtxy(nl,npatch_urban_tbd+3) * &
                (1.0_rk8 - urbinp%wtlunit_roof(nl,n))/3.0_rk8 * &
                (1.0_rk8 - urbinp%wtroad_perv(nl,n))
        wtxy(nl,npatch_urban_tbd+4) = wtxy(nl,npatch_urban_tbd+4) * &
                (1.0_rk8 - urbinp%wtlunit_roof(nl,n))/3.0_rk8 * &
                 urbinp%wtroad_perv(nl,n)
      end if

      ! Initialize urban high density weights
      n = udens_hd - udens_base
      do nurb = npatch_urban_hd , npatch_urban_md-1
        vegxy(nl,nurb) = noveg
        wtxy(nl,nurb)  = pcturb(nl,n) / 100._rk8
      end do
      if ( pcturb(nl,n) > 0.0_rk8 ) then
        wtxy(nl,npatch_urban_hd) = wtxy(nl,npatch_urban_hd) * &
                urbinp%wtlunit_roof(nl,n)
        wtxy(nl,npatch_urban_hd+1) = wtxy(nl,npatch_urban_hd+1) * &
                (1.0_rk8 - urbinp%wtlunit_roof(nl,n))/3.0_rk8
        wtxy(nl,npatch_urban_hd+2) = wtxy(nl,npatch_urban_hd+2) * &
                (1.0_rk8 - urbinp%wtlunit_roof(nl,n))/3.0_rk8
        wtxy(nl,npatch_urban_hd+3) = wtxy(nl,npatch_urban_hd+3) * &
                (1.0_rk8 - urbinp%wtlunit_roof(nl,n))/3.0_rk8 * &
                (1.0_rk8 - urbinp%wtroad_perv(nl,n))
        wtxy(nl,npatch_urban_hd+4) = wtxy(nl,npatch_urban_hd+4) * &
                (1.0_rk8 - urbinp%wtlunit_roof(nl,n))/3.0_rk8 * &
                 urbinp%wtroad_perv(nl,n)
      end if

      ! Initialize urban medium density weights
      n = udens_md - udens_base
      do nurb = npatch_urban_md , npatch_lake-1
        vegxy(nl,nurb) = noveg
        wtxy(nl,nurb)  = pcturb(nl,n) / 100._rk8
      end do
      if ( pcturb(nl,n) > 0.0_rk8 )then
        wtxy(nl,npatch_urban_md) = wtxy(nl,npatch_urban_md) * &
                urbinp%wtlunit_roof(nl,n)
        wtxy(nl,npatch_urban_md+1) = wtxy(nl,npatch_urban_md+1) * &
                (1.0_rk8 - urbinp%wtlunit_roof(nl,n))/3.0_rk8
        wtxy(nl,npatch_urban_md+2) = wtxy(nl,npatch_urban_md+2) * &
                (1.0_rk8 - urbinp%wtlunit_roof(nl,n))/3.0_rk8
        wtxy(nl,npatch_urban_md+3) = wtxy(nl,npatch_urban_md+3) * &
                (1.0_rk8 - urbinp%wtlunit_roof(nl,n))/3.0_rk8 * &
                (1.0_rk8 - urbinp%wtroad_perv(nl,n))
        wtxy(nl,npatch_urban_md+4) = wtxy(nl,npatch_urban_md+4) * &
                (1.0_rk8 - urbinp%wtlunit_roof(nl,n))/3.0_rk8 * &
                 urbinp%wtroad_perv(nl,n)
      end if
    end do

    ! Check to make sure we have valid urban data for each urban patch

    found = .false.
    do nl = begg , endg
      do n = 1 , numurbl
        if ( pcturb(nl,n) > 0.0_rk8 ) then
          if (urbinp%canyon_hwr(nl,n)            <= 0._rk8 .or. &
              urbinp%em_improad(nl,n)            <= 0._rk8 .or. &
              urbinp%em_perroad(nl,n)            <= 0._rk8 .or. &
              urbinp%em_roof(nl,n)               <= 0._rk8 .or. &
              urbinp%em_wall(nl,n)               <= 0._rk8 .or. &
              urbinp%ht_roof(nl,n)               <= 0._rk8 .or. &
              urbinp%thick_roof(nl,n)            <= 0._rk8 .or. &
              urbinp%thick_wall(nl,n)            <= 0._rk8 .or. &
              urbinp%t_building_max(nl,n)        <= 0._rk8 .or. &
              urbinp%t_building_min(nl,n)        <= 0._rk8 .or. &
              urbinp%wind_hgt_canyon(nl,n)       <= 0._rk8 .or. &
              urbinp%wtlunit_roof(nl,n)          <= 0._rk8 .or. &
              urbinp%wtroad_perv(nl,n)           <= 0._rk8 .or. &
              any(urbinp%alb_improad_dir(nl,n,:) <= 0._rk8) .or. &
              any(urbinp%alb_improad_dif(nl,n,:) <= 0._rk8) .or. &
              any(urbinp%alb_perroad_dir(nl,n,:) <= 0._rk8) .or. &
              any(urbinp%alb_perroad_dif(nl,n,:) <= 0._rk8) .or. &
              any(urbinp%alb_roof_dir(nl,n,:)    <= 0._rk8) .or. &
              any(urbinp%alb_roof_dif(nl,n,:)    <= 0._rk8) .or. &
              any(urbinp%alb_wall_dir(nl,n,:)    <= 0._rk8) .or. &
              any(urbinp%alb_wall_dif(nl,n,:)    <= 0._rk8) .or. &
              any(urbinp%tk_roof(nl,n,:)         <= 0._rk8) .or. &
              any(urbinp%tk_wall(nl,n,:)         <= 0._rk8) .or. &
              any(urbinp%cv_roof(nl,n,:)         <= 0._rk8) .or. &
              any(urbinp%cv_wall(nl,n,:)         <= 0._rk8)) then
            found = .true.
            nindx = nl
            dindx = n
            exit
          else
            if (urbinp%nlev_improad(nl,n) > 0) then
               nlev = urbinp%nlev_improad(nl,n)
               if (any(urbinp%tk_improad(nl,n,1:nlev) <= 0._rk8) .or. &
                   any(urbinp%cv_improad(nl,n,1:nlev) <= 0._rk8)) then
                  found = .true.
                  nindx = nl
                  dindx = n
                  exit
               end if
            end if
          end if
          if (found) exit
        end if
      end do
    end do
    if ( found ) then
      write(stderr,*)'surfrd error: no valid urban data for nl=',nindx
      write(stderr,*)'density type:    ',dindx
      write(stderr,*)'canyon_hwr:      ',urbinp%canyon_hwr(nindx,dindx)
      write(stderr,*)'em_improad:      ',urbinp%em_improad(nindx,dindx)
      write(stderr,*)'em_perroad:      ',urbinp%em_perroad(nindx,dindx)
      write(stderr,*)'em_roof:         ',urbinp%em_roof(nindx,dindx)
      write(stderr,*)'em_wall:         ',urbinp%em_wall(nindx,dindx)
      write(stderr,*)'ht_roof:         ',urbinp%ht_roof(nindx,dindx)
      write(stderr,*)'thick_roof:      ',urbinp%thick_roof(nindx,dindx)
      write(stderr,*)'thick_wall:      ',urbinp%thick_wall(nindx,dindx)
      write(stderr,*)'t_building_max:  ',urbinp%t_building_max(nindx,dindx)
      write(stderr,*)'t_building_min:  ',urbinp%t_building_min(nindx,dindx)
      write(stderr,*)'wind_hgt_canyon: ',urbinp%wind_hgt_canyon(nindx,dindx)
      write(stderr,*)'wtlunit_roof:    ',urbinp%wtlunit_roof(nindx,dindx)
      write(stderr,*)'wtroad_perv:     ',urbinp%wtroad_perv(nindx,dindx)
      write(stderr,*)'alb_improad_dir: ',urbinp%alb_improad_dir(nindx,dindx,:)
      write(stderr,*)'alb_improad_dif: ',urbinp%alb_improad_dif(nindx,dindx,:)
      write(stderr,*)'alb_perroad_dir: ',urbinp%alb_perroad_dir(nindx,dindx,:)
      write(stderr,*)'alb_perroad_dif: ',urbinp%alb_perroad_dif(nindx,dindx,:)
      write(stderr,*)'alb_roof_dir:    ',urbinp%alb_roof_dir(nindx,dindx,:)
      write(stderr,*)'alb_roof_dif:    ',urbinp%alb_roof_dif(nindx,dindx,:)
      write(stderr,*)'alb_wall_dir:    ',urbinp%alb_wall_dir(nindx,dindx,:)
      write(stderr,*)'alb_wall_dif:    ',urbinp%alb_wall_dif(nindx,dindx,:)
      write(stderr,*)'tk_roof:         ',urbinp%tk_roof(nindx,dindx,:)
      write(stderr,*)'tk_wall:         ',urbinp%tk_wall(nindx,dindx,:)
      write(stderr,*)'cv_roof:         ',urbinp%cv_roof(nindx,dindx,:)
      write(stderr,*)'cv_wall:         ',urbinp%cv_wall(nindx,dindx,:)
      if (urbinp%nlev_improad(nindx,dindx) > 0) then
        nlev = urbinp%nlev_improad(nindx,dindx)
        write(stderr,*)'tk_improad: ',urbinp%tk_improad(nindx,dindx,1:nlev)
        write(stderr,*)'cv_improad: ',urbinp%cv_improad(nindx,dindx,1:nlev)
      end if
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if
    do nl = begg , endg
      if ( pctspec(nl)-100._rk8 >  1.0e-4_rk8 ) then
        ldomain%pftm(nl) = 0
      end if
    end do
    deallocate(pctgla,pctlak,pctwet,pcturb,pcturb_tot)
  end subroutine surfrd_wtxy_special
  !
  ! Determine wtxy and veg arrays for non-dynamic landuse mode
  !
  subroutine surfrd_wtxy_veg_all(ncid,ldomain)
    use mod_clm_varctl , only : create_crop_landunit , irrigate
    use mod_clm_pftvarcon , only : nc3crop , nc3irrig , npcropmin ,  &
          ncorn , ncornirrig , nsoybean , nsoybeanirrig , nscereal , &
          nscerealirrig , nwcereal , nwcerealirrig , noveg
    use mod_clm_domain , only : domain_type
    implicit none
    type(clm_filetype) , intent(inout) :: ncid   ! netcdf id
    type(domain_type) , intent(inout) :: ldomain
    integer(ik4) :: m , nl         ! indices
    integer(ik4) :: begg , endg    ! beg/end gcell index
    integer(ik4) :: ier , tot_ier  ! error status
    real(rk8) :: sumpct            ! sum of %pft over maxpatch_pft
    ! percent of vegetated gridcell area for PFTs
    real(rk8) , allocatable , dimension(:,:) :: pctpft
    logical  :: crop = .false.  ! if crop data on this section of file
    character(len=32) :: subname = 'surfrd_wtxy_veg_all'  ! subroutine name

    call get_proc_bounds(begg,endg)
    allocate(pctpft(begg:endg,0:numpft))

    if ( .not. clm_check_dimlen(ncid, 'lsmpft', numpft+1) ) then
      write (stderr,*) 'LSMPFT in file not ', numpft+1
      call fatal(__FILE__,__LINE__, trim(subname)//' lsmpft /= numpft+1' )
    end if

    if ( clm_check_var(ncid,'PCT_PFT') ) then
      call clm_readvar(ncid,'PCT_PFT',pctpft(:,0:numpft),gcomm_gridcell)
    else
      call fatal(__FILE__,__LINE__, &
        trim(subname)//' ERROR: PCT_PFT NOT on surfdata file' )
    end if

    where (pctpft < eps_fact*epsilon(1.0) )
      pctpft = 0.0_rk8
    end where

    ier = 0
    do nl = begg , endg
      if ( ldomain%pftm(nl) >= 0 ) then
        ! Error check: make sure PFTs sum to 100% cover for vegetated landunit
        ! (convert pctpft from percent with respect to gridcel to percent with
        ! respect to vegetated landunit)
        ! THESE CHECKS NEEDS TO BE THE SAME AS IN pftdynMod.F90!
        if ( pctspec(nl)- 100._rk8 > 1.0e-4_rk8 ) then
          ! pctspec not within eps_fact*epsilon of 100
          sumpct = 0._rk8
          do m = 0 , numpft
            sumpct = sumpct + pctpft(nl,m) * 100._rk8/(100._rk8-pctspec(nl))
          end do
          if ( abs(sumpct - 100._rk8) > 1.0e-4_rk8 ) then
            write(stderr,*) 'SUMPFT  = ',sum(pctpft(nl,:))
            write(stderr,*) 'PCTSPEC = ',pctspec(nl)
            write(stderr,*) 'SUMPCT  = ',sumpct
            write(stderr,*) trim(subname)// &
                    ' ERROR: sum(pct) over numpft+1 is not = 100.'
            write(stderr,*) 'Unbalance ', sumpct-100._rk8, ' AT LAT,LON ', &
                    ldomain%latc(nl), ldomain%lonc(nl)
            ier = ier + 1
          end if
          if ( sumpct < -0.000001_rk8 ) then
            write(stderr,*) trim(subname)// &
                    ' ERROR: sum(pct) over numpft+1 is < 0.'
            write(stderr,*) 'SUMPCT = ', sumpct, ' at nl : ',nl
            ier = ier + 100000
          end if
          do m = npcropmin , numpft
#ifdef CNDV
            if ( pctpft(nl,m) > 0.0_rk8 ) then
              pctpft(nl,noveg) = pctpft(nl,noveg) + pctpft(nl,m)
              pctpft(nl,m) = 0.0_rk8
            end if
#else
            if ( pctpft(nl,m) > 0.0_rk8 ) crop = .true.
#endif
          end do
        end if
        ! Set weight of each pft wrt gridcell
        ! (note that maxpatch_pft = numpft+1 here)
        do m = 1 , numpft+1
          vegxy(nl,m) = m - 1 ! 0 (bare ground) to numpft
          wtxy(nl,m) = pctpft(nl,m-1) / 100._rk8
        end do
      end if
    end do

#ifdef CROP
    crop = .true.
#endif
    call sumall(ier,tot_ier)
    if ( tot_ier /= 0 ) then
      write(stderr,*) 'Total errors = ', tot_ier
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    call trueforall(crop,crop_prog)
    if ( crop_prog .and. .not. create_crop_landunit ) then
      call fatal(__FILE__,__LINE__, &
              trim(subname)//' ERROR: prognostic crop '// &
              'PFTs require create_crop_landunit=.true.' )
    end if
#ifdef DYNPFT
    if ( crop_prog ) then
      call fatal(__FILE__,__LINE__, &
              trim(subname)//' ERROR: prognostic crop '// &
              'is incompatible with transient landuse' )
    end if
#endif
    if ( .not. crop_prog .and. irrigate ) then
      call fatal(__FILE__,__LINE__,trim(subname)// &
           ' ERROR surfrdMod: irrigate = .true. requires CROP model active.' )
    end if

    if ( myid == italk .and. crop_prog .and. .not. irrigate ) then
       write(stdout,*) trim(subname)// &
         ' crop=.T. and irrigate=.F., so merging irrigated pfts with rainfed'
         ! in the following do-loop
    end if

    ! repeat do-loop for error checking and for rainfed crop case

    do nl = begg , endg
      if ( ldomain%pftm(nl) >= 0 ) then
        if ( pctspec(nl) - 100._rk8 >  1.0e-4_rk8 ) then
          ! pctspec not within eps_fact*epsilon of 100
          if ( .not. crop_prog .and. &
                wtxy(nl,nc3irrig+1) > eps_fact*epsilon(1.0) ) then
            write (stderr,*) 'wtxy(nl,nc3irrig+1) = ', wtxy(nl,nc3irrig+1)
            call fatal(__FILE__,__LINE__, &
                trim(subname)//' ERROR surfrdMod: irrigated crop PFT '//&
                'requires CROP model active.' )
          end if
          if ( crop_prog .and. .not. irrigate ) then
            wtxy(nl,nc3crop+1) = wtxy(nl,nc3crop+1) + wtxy(nl,nc3irrig+1)
            wtxy(nl,nc3irrig+1) = 0._rk8
            wtxy(nl,ncorn+1) = wtxy(nl,ncorn+1) + wtxy(nl,ncornirrig+1)
            wtxy(nl,ncornirrig+1) = 0._rk8
            wtxy(nl,nscereal+1) = wtxy(nl,nscereal+1) + wtxy(nl,nscerealirrig+1)
            wtxy(nl,nscerealirrig+1) = 0._rk8
            wtxy(nl,nwcereal+1) = wtxy(nl,nwcereal+1) + wtxy(nl,nwcerealirrig+1)
            wtxy(nl,nwcerealirrig+1) = 0._rk8
            wtxy(nl,nsoybean+1) = wtxy(nl,nsoybean+1) + wtxy(nl,nsoybeanirrig+1)
            wtxy(nl,nsoybeanirrig+1) = 0._rk8
          end if
        end if
      end if
    end do
    deallocate(pctpft)
  end subroutine surfrd_wtxy_veg_all
  !
  ! Determine wtxy and vegxy for CNDV mode.
  !
  subroutine surfrd_wtxy_veg_dgvm()
    use mod_clm_pftvarcon , only : noveg , crop
    use mod_clm_varctl , only : create_crop_landunit
    implicit none
    integer(ik4) :: m , nl        ! indices
    integer(ik4) :: begg , endg   ! beg/end gcell index

    call get_proc_bounds(begg,endg)
    do nl = begg , endg
      do m = 1 , maxpatch_pft ! CNDV means allocate_all_vegpfts = .true.
        if ( create_crop_landunit ) then ! been through surfrd_wtxy_veg_all
          if ( crop(m-1) == 0 ) then     ! so update natural vegetation only
            wtxy(nl,m) = 0._rk8          ! crops should have values >= 0.
          end if
        else                    ! not been through surfrd_wtxy_veg_all
          wtxy(nl,m) = 0._rk8   ! so update all vegetation
          vegxy(nl,m) = m - 1 ! 0 (bare ground) to maxpatch_pft-1 (= 16)
        end if
      end do
      ! bare ground weights
      wtxy(nl,noveg+1) = max(0._rk8, 1._rk8 - sum(wtxy(nl,:)))
      if (abs(sum(wtxy(nl,:)) - 1._rk8) > 1.e-5_rk8) then
        write(stderr,*) 'all wtxy =', wtxy(nl,:)
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if ! error check
    end do
  end subroutine surfrd_wtxy_veg_dgvm

end module mod_clm_surfrd
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
