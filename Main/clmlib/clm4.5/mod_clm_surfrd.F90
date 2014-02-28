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
  use mod_dynparam , only : myid , ds
  use mod_mppparam
  use mod_clm_nchelper
  use mod_clm_varpar , only : nlevsoifl , numpft , maxpatch_pft , numcft , &
         maxpatch , npatch_urban_tbd , npatch_urban_hd , npatch_urban_md , &
         numurbl , npatch_lake , npatch_wet , npatch_glacier ,             &
         maxpatch_urb , npatch_glacier_mec , maxpatch_glcmec
  use mod_clm_varctl , only : glc_topomax
  use mod_clm_varsur , only : wtxy , vegxy , pctspec
  use mod_clm_decomp , only : get_proc_bounds , procinfo , numg
  use mod_clm_type

  implicit none

  private

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
  real(rk8) , private, parameter :: eps_fact = 2.D0

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
    implicit none
    type(domain_type) , intent(inout) :: ldomain   ! domain to init
    character(len=*) , intent(in)    :: filename  ! grid filename
    type(clm_filetype) :: ncid       ! netcdf id
    integer(ik4) :: ibeg             ! local beg index
    integer(ik4) :: iend             ! local end index
    integer(ik4) :: ni , nj , ns     ! size of grid on file
    integer(ik4) :: inni , innj      ! size of grid on file
    integer(ik4) :: ier , ret        ! error status
    real(rk8) , allocatable , dimension(:,:) :: rdata2d ! temporary
    real(rk8) , allocatable , dimension(:) :: rdata1d ! temporary
    integer(ik4) :: n , i , j        ! indices
    real(rk8) :: eps = 1.0D-12       ! lat/lon error tolerance
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

    ! RegCM INTERNAL grid is 2:jx-2,2:iy-2
    allocate(rdata2d(inni,innj))
    allocate(rdata1d(numg))
    ni = inni - 3
    nj = innj - 3
    ns = ni*nj

    call get_proc_bounds(ibeg,iend)
    call domain_init(ldomain,ni=ni,nj=nj,nbeg=ibeg,nend=iend)

    ! Read mask from input file
    call clm_readvar(ncid,'mask',rdata2d)
    call getmem2d(procinfo%gcmask,1,ni,1,nj,'surfrd:gcmask')
    procinfo%gcmask = .false.
    do j = 2 , innj-2
      do i = 2 , inni-2
        procinfo%gcmask(i-1,j-1) = (rdata2d(i,j) > 0)
      end do
    end do
    if ( count(procinfo%gcmask) /= numg ) then
      write(stderr,*) 'Unmatch from landmask in ATM and CLM'
      write(stderr,*) 'ATM : ', count(procinfo%gcmask)
      write(stderr,*) 'CLM : ', numg
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    rdata1d = pack(rdata2d(2:inni-2,2:innj-2),procinfo%gcmask)
    ldomain%frac = rdata1d(ibeg:iend)

    ! Read in lon, lat from RegCM Subgrid domain file
    call clm_readvar(ncid,'xlon',rdata2d)
    rdata1d = pack(rdata2d(2:inni-2,2:innj-2),procinfo%gcmask)
    ldomain%lonc = rdata1d(ibeg:iend)
    call clm_readvar(ncid,'xlat',rdata2d)
    rdata1d = pack(rdata2d(2:inni-2,2:innj-2),procinfo%gcmask)
    ldomain%latc = rdata1d(ibeg:iend)
    call clm_readvar(ncid,'topo',rdata2d)
    rdata1d = pack(rdata2d(2:inni-2,2:innj-2),procinfo%gcmask)
    ldomain%topo = rdata1d(ibeg:iend)

    deallocate(rdata2d)
    deallocate(rdata1d)
    call clm_closefile(ncid)

    where ( ldomain%frac > 0 ) ldomain%frac = 100.0D0
    ldomain%area = ds*ds
    where ( ldomain%frac > 0 )
      ldomain%mask = 1
    else where
      ldomain%mask = 0
    end where
    ldomain%pftm = ldomain%mask

    ! Check lat limited to -90,90

    if (minval(ldomain%latc) < -90.0D0 .or. &
        maxval(ldomain%latc) >  90.0D0) then
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
    use mod_clm_varctl , only : allocate_all_vegpfts , create_crop_landunit
    use mod_clm_pftvarcon , only : noveg
    use mod_clm_domain , only : domain_type , domain_init , domain_clean
    implicit none
    ! land domain associated with wtxy
    type(domain_type) , intent(inout) :: ldomain
    character(len=*) , intent(inout) :: lfsurdat ! surface dataset filename
    ! local domain associated with surface dataset
    type(domain_type) :: surfdata_domain
    integer(ik4) :: n                        ! loop indices
    integer(ik4) :: ni , nj , ns             ! domain sizes
    integer(ik4) :: inni , innj              ! domain sizes
    logical :: readvar                       ! true => variable is on dataset
    real(rk8) :: rmaxlon , rmaxlat           ! local min/max vars
    type(clm_filetype) :: ncid               ! netcdf id
    integer(ik4) :: begg , endg              ! beg,end gridcell indices
    logical :: istype_domain              ! true => input file is of type domain
    real(rk8) , allocatable , dimension(:,:) :: rdata2d ! temporary
    real(rk8) , allocatable , dimension(:) :: rdata1d ! temporary
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

    vegxy(:,:) = noveg
    wtxy(:,:)  = 0.D0
    pctspec(:) = 0.D0

    ! Read surface data

    call clm_openfile(lfsurdat,ncid)

    ! Determine dimensions

    call clm_inqdim(ncid,'jx',inni)
    call clm_inqdim(ncid,'iy',innj)

    ! RegCM INTERNAL grid is 2:jx-2,2:iy-2
    allocate(rdata2d(inni,innj))
    allocate(rdata1d(numg))
    ni = inni - 3
    nj = innj - 3
    ns = ni*nj

    call domain_init(surfdata_domain,ni,nj,begg,endg,clmlevel=grlnd)

    if ( ldomain%ni /= surfdata_domain%ni .or. &
         ldomain%nj /= surfdata_domain%nj ) then
      write(stderr,*) &
              trim(subname)//': surfdata/lon/lat size mismatch error'
      call fatal(__FILE__,__LINE__,trim(subname))
    end if
    call domain_clean(surfdata_domain)

    ! Obtain special landunit info

    call surfrd_wtxy_special(ncid, inni,innj,rdata2d,rdata1d)

    ! Obtain vegetated landunit info

#if (defined CNDV)
    if (create_crop_landunit) then ! CNDV means allocate_all_vegpfts = .true.
      call surfrd_wtxy_veg_all(ncid,inni,innj,rdata2d,rdata1d,ldomain)
    end if
    call surfrd_wtxy_veg_dgvm()
#else
    if (allocate_all_vegpfts) then
      call surfrd_wtxy_veg_all(ncid,inni,innj,rdata2d,rdata1d,ldomain)
    else
      call fatal(__FILE__,__LINE__, &
              trim(subname) // 'only allocate_all_vegpfts is supported')
    end if
#endif
    call clm_closefile(ncid)
    deallocate(rdata2d,rdata1d)

    if ( myid == italk )then
      write(stdout,*) 'Successfully read surface boundary data'
      write(stdout,*)
    end if

  end subroutine surfrd_get_data
  !
  ! Determine weight with respect to gridcell of all special "pfts" as well
  ! as soil color and percent sand and clay
  !
  subroutine surfrd_wtxy_special(ncid,inni,innj,rdata2d,rdata1d)
    use mod_clm_pftvarcon , only : noveg
    use mod_clm_urbaninput , only : urbinp
    use mod_clm_varpar , only : maxpatch_glcmec , nlevurb
    use mod_clm_varcon , only : udens_base , udens_tbd , udens_hd , udens_md 
    implicit none
    type(clm_filetype) , intent(inout) :: ncid  ! netcdf id
    real(rk8) , intent(inout) , dimension(:,:) :: rdata2d
    real(rk8) , intent(inout) , dimension(:) :: rdata1d
    integer(ik4) , intent(in) :: inni , innj
    integer(ik4)  :: n , nl , nurb , g          ! indices
    integer(ik4)  :: begg , endg                ! gcell beg/end
    integer(ik4)  :: dimid , varid              ! netCDF id's
    real(rk8) , dimension(nlevsoifl) :: nlevsoidata
    logical :: found                      ! temporary for error check
    integer(ik4) :: nindx                 ! temporary for error check
    integer(ik4) :: ier                   ! error status
    integer(ik4) :: nlev                  ! level
    integer(ik4) :: dindx                 ! temporary for error check
    integer(ik4) :: npatch
    logical :: readvar
    real(rk8) :: closelat , closelon
    ! percent of grid cell is glacier
    real(rk8) , pointer , dimension(:) :: pctgla
    ! percent of grid cell is lake
    real(rk8) , pointer , dimension(:) :: pctlak
    ! percent of grid cell is wetland
    real(rk8) , pointer , dimension(:) :: pctwet
    ! percent of grid cell is urbanized
    real(rk8) , pointer , dimension(:,:) :: pcturb
    ! percent of grid cell is glacier_mec (in each elev class)
    real(rk8) , pointer , dimension(:,:) :: pctglc_mec
    ! percent of grid cell is glacier (sum over classes)
    real(rk8) , pointer , dimension(:) :: pctglc_mec_tot
    ! percent of grid cell is urban (sum over density classes)
    real(rk8) , pointer , dimension(:) :: pcturb_tot
    ! surface elevation in each elev class
    real(rk8) , pointer , dimension(:,:) :: topoglc_mec
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
      call clm_readvar(ncid,'PCT_WETLAND',rdata2d)
      rdata1d = pack(rdata2d(2:inni-2,2:innj-2),procinfo%gcmask)
      pctwet = rdata1d(begg:endg)
    else
      call fatal(__FILE__,__LINE__, &
        trim(subname)//' ERROR: PCT_WETLAND  NOT on surfdata file' )
    end if

    if ( clm_check_var(ncid,'PCT_LAKE') ) then
      call clm_readvar(ncid,'PCT_LAKE',rdata2d)
      rdata1d = pack(rdata2d(2:inni-2,2:innj-2),procinfo%gcmask)
      pctlak = rdata1d(begg:endg)
    else
      call fatal(__FILE__,__LINE__, &
        trim(subname)//' ERROR: PCT_LAKE NOT on surfdata file' )
    end if

    if ( clm_check_var(ncid,'PCT_GLACIER') ) then
      call clm_readvar(ncid,'PCT_GLACIER',rdata2d)
      rdata1d = pack(rdata2d(2:inni-2,2:innj-2),procinfo%gcmask)
      pctgla = rdata1d(begg:endg)
    else
      call fatal(__FILE__,__LINE__, &
        trim(subname)//' ERROR: PCT_GLACIER NOT on surfdata file' )
    end if

    ! If PCT_URBAN is not multi-density then set pcturb and nlevurb to zero 
    if (nlevurb == 0) then
      pcturb = 0.D0
      if ( myid == italk ) &
        write(stdout,*)'PCT_URBAN is not multi-density, pcturb set to 0'
    else
      if ( clm_check_var(ncid,'PCT_URBAN') ) then
        do n = 1 , numurbl
          call clm_readvar(ncid,'PCT_URBAN',rdata2d,n)
          rdata1d = pack(rdata2d(2:inni-2,2:innj-2),procinfo%gcmask)
          pcturb(n,:) = rdata1d(begg:endg)
        end do
      else
        call fatal(__FILE__,__LINE__, &
          trim(subname)//' ERROR: PCT_URBAN NOT on surfdata file' )
      end if
    end if
    if ( nlevurb == 0 )then
      if ( any(pcturb > 0.0D0) ) then
        call fatal(__FILE__,__LINE__, &
          trim(subname)//' ERROR: PCT_URBAN MUST be zero when nlevurb=0' )
      end if
    end if

    pcturb_tot(:) = 0.D0
    do n = 1 , numurbl
      do nl = begg , endg
        pcturb_tot(nl) = pcturb_tot(nl) + pcturb(nl,n)
      end do
    end do

    pctspec = pctwet + pctlak + pcturb_tot + pctgla

    ! Error check: glacier, lake, wetland, urban sum must be less than 100

    found = .false.
    do nl = begg , endg
      if (pctspec(nl) > 100.D0+1.D-04) then
        found = .true.
        nindx = nl
        exit
      end if
      if (found) exit
    end do
    if ( found ) then
      write(stderr,*)'surfrd error: PFT cover>100 for nl=',nindx
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! Determine veg and wtxy for special landunits

    do nl = begg , endg
      vegxy(nl,npatch_lake)   = noveg
      wtxy(nl,npatch_lake)    = pctlak(nl)/100.D0

      vegxy(nl,npatch_wet)    = noveg
      wtxy(nl,npatch_wet)     = pctwet(nl)/100.D0

      vegxy(nl,npatch_glacier)= noveg
      wtxy(nl,npatch_glacier) = pctgla(nl)/100.D0

      ! Initialize urban tall building district weights
      n = udens_tbd - udens_base
      do nurb = npatch_urban_tbd , npatch_urban_hd-1 
        vegxy(nl,nurb) = noveg
        wtxy(nl,nurb)  = pcturb(nl,n) / 100.D0
      end do
      if ( pcturb(nl,n) > 0.0D0 )then
        wtxy(nl,npatch_urban_tbd) = wtxy(nl,npatch_urban_tbd) * &
                urbinp%wtlunit_roof(nl,n)
        wtxy(nl,npatch_urban_tbd+1) = wtxy(nl,npatch_urban_tbd+1) * &
                (1.0D0 - urbinp%wtlunit_roof(nl,n))/3.0D0
        wtxy(nl,npatch_urban_tbd+2) = wtxy(nl,npatch_urban_tbd+2) * &
                (1.0D0 - urbinp%wtlunit_roof(nl,n))/3.0D0
        wtxy(nl,npatch_urban_tbd+3) = wtxy(nl,npatch_urban_tbd+3) * &
                (1.0D0 - urbinp%wtlunit_roof(nl,n))/3.0D0 * &
                (1.0D0 - urbinp%wtroad_perv(nl,n))
        wtxy(nl,npatch_urban_tbd+4) = wtxy(nl,npatch_urban_tbd+4) * &
                (1.0D0 - urbinp%wtlunit_roof(nl,n))/3.0D0 * &
                 urbinp%wtroad_perv(nl,n)
      end if

      ! Initialize urban high density weights
      n = udens_hd - udens_base
      do nurb = npatch_urban_hd , npatch_urban_md-1 
        vegxy(nl,nurb) = noveg
        wtxy(nl,nurb)  = pcturb(nl,n) / 100.D0
      end do
      if ( pcturb(nl,n) > 0.0D0 ) then
        wtxy(nl,npatch_urban_hd) = wtxy(nl,npatch_urban_hd) * &
                urbinp%wtlunit_roof(nl,n)
        wtxy(nl,npatch_urban_hd+1) = wtxy(nl,npatch_urban_hd+1) * &
                (1.0D0 - urbinp%wtlunit_roof(nl,n))/3.0D0
        wtxy(nl,npatch_urban_hd+2) = wtxy(nl,npatch_urban_hd+2) * &
                (1.0D0 - urbinp%wtlunit_roof(nl,n))/3.0D0
        wtxy(nl,npatch_urban_hd+3) = wtxy(nl,npatch_urban_hd+3) * &
                (1.0D0 - urbinp%wtlunit_roof(nl,n))/3.0D0 * &
                (1.0D0 - urbinp%wtroad_perv(nl,n))
        wtxy(nl,npatch_urban_hd+4) = wtxy(nl,npatch_urban_hd+4) * &
                (1.0D0 - urbinp%wtlunit_roof(nl,n))/3.0D0 * &
                 urbinp%wtroad_perv(nl,n)
      end if

      ! Initialize urban medium density weights
      n = udens_md - udens_base
      do nurb = npatch_urban_md , npatch_lake-1 
        vegxy(nl,nurb) = noveg
        wtxy(nl,nurb)  = pcturb(nl,n) / 100.D0
      end do
      if ( pcturb(nl,n) > 0.0D0 )then
        wtxy(nl,npatch_urban_md) = wtxy(nl,npatch_urban_md) * &
                urbinp%wtlunit_roof(nl,n)
        wtxy(nl,npatch_urban_md+1) = wtxy(nl,npatch_urban_md+1) * &
                (1.0D0 - urbinp%wtlunit_roof(nl,n))/3.0D0
        wtxy(nl,npatch_urban_md+2) = wtxy(nl,npatch_urban_md+2) * &
                (1.0D0 - urbinp%wtlunit_roof(nl,n))/3.0D0
        wtxy(nl,npatch_urban_md+3) = wtxy(nl,npatch_urban_md+3) * &
                (1.0D0 - urbinp%wtlunit_roof(nl,n))/3.0D0 * &
                (1.0D0 - urbinp%wtroad_perv(nl,n))
        wtxy(nl,npatch_urban_md+4) = wtxy(nl,npatch_urban_md+4) * &
                (1.0D0 - urbinp%wtlunit_roof(nl,n))/3.0D0 * &
                 urbinp%wtroad_perv(nl,n)
      end if
    end do

    ! Check to make sure we have valid urban data for each urban patch

    found = .false.
    do nl = begg , endg
      do n = 1 , numurbl
        if ( pcturb(nl,n) > 0.0D0 ) then
          if (urbinp%canyon_hwr(nl,n)            .le. 0.D0 .or. &
              urbinp%em_improad(nl,n)            .le. 0.D0 .or. &
              urbinp%em_perroad(nl,n)            .le. 0.D0 .or. &
              urbinp%em_roof(nl,n)               .le. 0.D0 .or. &
              urbinp%em_wall(nl,n)               .le. 0.D0 .or. &
              urbinp%ht_roof(nl,n)               .le. 0.D0 .or. &
              urbinp%thick_roof(nl,n)            .le. 0.D0 .or. &
              urbinp%thick_wall(nl,n)            .le. 0.D0 .or. &
              urbinp%t_building_max(nl,n)        .le. 0.D0 .or. &
              urbinp%t_building_min(nl,n)        .le. 0.D0 .or. &
              urbinp%wind_hgt_canyon(nl,n)       .le. 0.D0 .or. &
              urbinp%wtlunit_roof(nl,n)          .le. 0.D0 .or. &
              urbinp%wtroad_perv(nl,n)           .le. 0.D0 .or. &
              any(urbinp%alb_improad_dir(nl,n,:) .le. 0.D0) .or. &
              any(urbinp%alb_improad_dif(nl,n,:) .le. 0.D0) .or. &
              any(urbinp%alb_perroad_dir(nl,n,:) .le. 0.D0) .or. &
              any(urbinp%alb_perroad_dif(nl,n,:) .le. 0.D0) .or. &
              any(urbinp%alb_roof_dir(nl,n,:)    .le. 0.D0) .or. &
              any(urbinp%alb_roof_dif(nl,n,:)    .le. 0.D0) .or. &
              any(urbinp%alb_wall_dir(nl,n,:)    .le. 0.D0) .or. &
              any(urbinp%alb_wall_dif(nl,n,:)    .le. 0.D0) .or. &
              any(urbinp%tk_roof(nl,n,:)         .le. 0.D0) .or. &
              any(urbinp%tk_wall(nl,n,:)         .le. 0.D0) .or. &
              any(urbinp%cv_roof(nl,n,:)         .le. 0.D0) .or. &
              any(urbinp%cv_wall(nl,n,:)         .le. 0.D0)) then
            found = .true.
            nindx = nl
            dindx = n
            exit
          else
            if (urbinp%nlev_improad(nl,n) .gt. 0) then
               nlev = urbinp%nlev_improad(nl,n)
               if (any(urbinp%tk_improad(nl,n,1:nlev) .le. 0.D0) .or. &
                   any(urbinp%cv_improad(nl,n,1:nlev) .le. 0.D0)) then
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
      if (urbinp%nlev_improad(nindx,dindx) .gt. 0) then
        nlev = urbinp%nlev_improad(nindx,dindx)
        write(stderr,*)'tk_improad: ',urbinp%tk_improad(nindx,dindx,1:nlev)
        write(stderr,*)'cv_improad: ',urbinp%cv_improad(nindx,dindx,1:nlev)
      end if
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    deallocate(pctgla,pctlak,pctwet,pcturb,pcturb_tot)
  end subroutine surfrd_wtxy_special
  !
  ! Determine wtxy and veg arrays for non-dynamic landuse mode
  !
  subroutine surfrd_wtxy_veg_all(ncid,inni,innj,rdata2d,rdata1d,ldomain)
    use mod_clm_varctl , only : create_crop_landunit , fpftdyn , irrigate
    use mod_clm_pftvarcon , only : nc3crop , nc3irrig , npcropmin ,  &
          ncorn , ncornirrig , nsoybean , nsoybeanirrig , nscereal , &
          nscerealirrig , nwcereal , nwcerealirrig
    use mod_clm_domain , only : domain_type
    implicit none
    type(clm_filetype) , intent(inout) :: ncid   ! netcdf id
    real(rk8) , intent(inout) , dimension(:,:) :: rdata2d
    real(rk8) , intent(inout) , dimension(:) :: rdata1d
    integer(ik4) , intent(in) :: inni , innj
    type(domain_type) , intent(inout) :: ldomain
    integer(ik4) :: m , mp7 , mp8 , mp11 , n , nl ! indices
    integer(ik4) :: begg , endg               ! beg/end gcell index
    integer(ik4) :: dimid , varid             ! netCDF id's
    integer(ik4) :: ier                       ! error status  
    logical :: readvar                        ! is variable on dataset
    real(rk8) :: sumpct                       ! sum of %pft over maxpatch_pft
    ! percent of vegetated gridcell area for PFTs
    real(rk8) , allocatable , dimension(:,:) :: pctpft
    real(rk8) , dimension(0:numpft) :: numpftp1data
    logical  :: crop = .false.  ! if crop data on this section of file
    character(len=32) :: subname = 'surfrd_wtxy_veg_all'  ! subroutine name

    call get_proc_bounds(begg,endg)
    allocate(pctpft(begg:endg,0:numpft))

    if ( .not. clm_check_dimlen(ncid, 'lsmpft', numpft+1) ) then
      call fatal(__FILE__,__LINE__, trim(subname)//' lsmpft /= numpft+1' )
    end if

    if ( clm_check_var(ncid,'PCT_PFT') ) then
      do n = 1 , numpft
        call clm_readvar(ncid,'PCT_PFT',rdata2d,n)
        rdata1d = pack(rdata2d(2:inni-2,2:innj-2),procinfo%gcmask)
        pctpft(:,n) = rdata1d(begg:endg)
      end do
    else
      call fatal(__FILE__,__LINE__, &
        trim(subname)//' ERROR: PCT_PFT NOT on surfdata file' )
    end if

    where (pctpft < eps_fact*epsilon(1.D0) )
      pctpft = 0.0D0
    end where

    do nl = begg , endg
      if ( ldomain%pftm(nl) >= 0 ) then
        ! Error check: make sure PFTs sum to 100% cover for vegetated landunit 
        ! (convert pctpft from percent with respect to gridcel to percent with 
        ! respect to vegetated landunit)
        ! THESE CHECKS NEEDS TO BE THE SAME AS IN pftdynMod.F90!
        if ( pctspec(nl) < 100.D0 * (1.D0 - eps_fact*epsilon(1.D0)) ) then
          ! pctspec not within eps_fact*epsilon of 100
          sumpct = 0.D0
          do m = 0 , numpft
            sumpct = sumpct + pctpft(nl,m) * 100.D0/(100.D0-pctspec(nl))
          end do
          if ( abs(sumpct - 100.D0) > 0.1D-4 ) then
            write(stderr,*) trim(subname)// &
                    ' ERROR: sum(pct) over numpft+1 is not = 100.'
            write(stderr,*) 'Unbalance ', sumpct-100.D0, ' AT LAT,LON ', &
                    ldomain%latc(nl), ldomain%lonc(nl)
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
          if ( sumpct < -0.000001D0 ) then
            write(stderr,*) trim(subname)// &
                    ' ERROR: sum(pct) over numpft+1 is < 0.'
            write(stderr,*) sumpct, nl
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
          do m = npcropmin , numpft
            if ( pctpft(nl,m) > 0.0D0 ) crop = .true.
          end do
        end if
        ! Set weight of each pft wrt gridcell
        ! (note that maxpatch_pft = numpft+1 here)
        do m = 1 , numpft+1
          vegxy(nl,m)  = m - 1 ! 0 (bare ground) to numpft
          wtxy(nl,m) = pctpft(nl,m-1) / 100.D0
        end do
      end if
    end do

    call trueforall(crop,crop_prog)
    if ( crop_prog .and. .not. create_crop_landunit ) then
      call fatal(__FILE__,__LINE__, &
              trim(subname)//' ERROR: prognostic crop '// &
              'PFTs require create_crop_landunit=.true.' )
    end if
    if ( crop_prog .and. fpftdyn /= ' ' ) then
      call fatal(__FILE__,__LINE__, &
              trim(subname)//' ERROR: prognostic crop '// &
              'is incompatible with transient landuse' )
    end if
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
        if ( pctspec(nl) < 100.D0 * (1.D0 - eps_fact*epsilon(1.D0)) ) then
          ! pctspec not within eps_fact*epsilon of 100
          if ( .not. crop_prog .and. wtxy(nl,nc3irrig+1) > 0.D0 ) then
            call fatal(__FILE__,__LINE__, &
                trim(subname)//' ERROR surfrdMod: irrigated crop PFT '//&
                'requires CROP model active.' )
          end if
          if ( crop_prog .and. .not. irrigate ) then
            wtxy(nl,nc3crop+1) = wtxy(nl,nc3crop+1) + wtxy(nl,nc3irrig+1)
            wtxy(nl,nc3irrig+1) = 0.D0
            wtxy(nl,ncorn+1) = wtxy(nl,ncorn+1) + wtxy(nl,ncornirrig+1)
            wtxy(nl,ncornirrig+1) = 0.D0
            wtxy(nl,nscereal+1) = wtxy(nl,nscereal+1) + wtxy(nl,nscerealirrig+1)
            wtxy(nl,nscerealirrig+1) = 0.D0
            wtxy(nl,nwcereal+1) = wtxy(nl,nwcereal+1) + wtxy(nl,nwcerealirrig+1)
            wtxy(nl,nwcerealirrig+1) = 0.D0
            wtxy(nl,nsoybean+1) = wtxy(nl,nsoybean+1) + wtxy(nl,nsoybeanirrig+1)
            wtxy(nl,nsoybeanirrig+1) = 0.D0
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
            wtxy(nl,m) = 0.D0       ! crops should have values >= 0.
          end if
        else                        ! not been through surfrd_wtxy_veg_all
          wtxy(nl,m) = 0.D0   ! so update all vegetation
          vegxy(nl,m) = m - 1 ! 0 (bare ground) to maxpatch_pft-1 (= 16)
        end if
      end do
      ! bare ground weights
      wtxy(nl,noveg+1) = max(0.D0, 1.D0 - sum(wtxy(nl,:)))
      if (abs(sum(wtxy(nl,:)) - 1.D0) > 1.D-5) then
        write(stderr,*) 'all wtxy =', wtxy(nl,:)
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if ! error check
    end do
  end subroutine surfrd_wtxy_veg_dgvm
   
end module mod_clm_surfrd
