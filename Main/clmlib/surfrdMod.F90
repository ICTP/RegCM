#include <misc.h>
#include <preproc.h>

module surfrdMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: surfrdMod
!
! !DESCRIPTION:
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
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_varpar  , only : lsmlon, lsmlat
  use clm_varpar  , only : nlevsoi, numpft, &
                           maxpatch_pft, maxpatch_cft, maxpatch, &
                           npatch_urban, npatch_lake, npatch_wet, npatch_glacier
  use clm_varsur  , only : wtxy, vegxy
  use clm_varsur  , only : pctspec , pctspecB
  use ncdio
  use clmtype
  use spmdMod                         
  use clm_varctl,   only : scmlat, scmlon, single_column
  use decompMod   , only : get_proc_bounds,gsMap_lnd_gdc2glo,perm_lnd_gdc2glo
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: surfrd  ! Read surface dataset and determine subgrid weights
  public :: surfrd_get_grid  ! Read surface dataset into domain
  public :: surfrd_get_latlon  ! Read surface dataset into domain
  public :: surfrd_get_frac  ! Read land fraction into domain
  public :: surfrd_get_topo  ! Read topography into domain
! abt rcm below
  public :: rcmsurfrd_get_grid     ! Read surface dataset into domain
  public :: rcmsurfrd_get_latlon   ! Read surface dataset into domain
  public :: rcmsurfrd_get_frac     ! Read land fraction into domain
  public :: rcmsurfrd_get_topo     ! Read topography into domain
  public :: rcmsurfrd_bvocs        ! Read biogenic emissions
  private :: pft_adjustment        ! If there are gridcells with no
                                   ! land cover type  
  private :: clm2bats_conversion   ! convert CLM types to BATS types
  public  :: clm_getsoitex         ! convert CLM soil types to BATS types
! abt rcm above
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Updated by T Craig
! Modified by Ahmed Tawfik for RegCM
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: surfrd_wtxy_special
  private :: surfrd_wtxy_veg_rank
  private :: surfrd_wtxy_veg_all
  private :: surfrd_wtxy_veg_dgvm
  private :: surfrd_mkrank
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd
!
! !INTERFACE:
!  subroutine surfrd(vegxy, wtxy, lfsurdat, domain)
  subroutine surfrd(lfsurdat, domain)
!
! !DESCRIPTION:
! Read the surface dataset and create subgrid weights.
! The model's surface dataset recognizes 5 basic land cover types within
! a grid cell: lake, wetland, urban, glacier, and vegetated. The vegetated
! portion of the grid cell is comprised of up to [maxpatch_pft] PFTs. These
! subgrid patches are read in explicitly for each grid cell. This is in
! contrast to LSMv1, where the PFTs were built implicitly from biome types.
!    o real edges of grid
!    o integer  number of longitudes per latitude
!    o real latitude  of grid cell (degrees)
!    o real longitude of grid cell (degrees)
!    o integer surface type: 0 = ocean or 1 = land
!    o integer soil color (1 to 20) for use with soil albedos
!    o real soil texture, %sand, for thermal and hydraulic properties
!    o real soil texture, %clay, for thermal and hydraulic properties
!    o real % of cell covered by lake    for use as subgrid patch
!    o real % of cell covered by wetland for use as subgrid patch
!    o real % of cell that is urban      for use as subgrid patch
!    o real % of cell that is glacier    for use as subgrid patch
!    o integer PFTs
!    o real % abundance PFTs (as a percent of vegetated area)
!
! !USES:
    use clm_varctl  , only : allocate_all_vegpfts
    use pftvarcon   , only : noveg
    use fileutils   , only : getfil
    use domainMod , only : domain_type
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(in) :: lfsurdat               ! surf filename
    type(domain_type),intent(in) :: domain ! domain associated with wtxy
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: locfn                          ! local file name
    integer  :: ncid,dimid,varid                         ! netCDF id's
    integer  :: begg,endg   
    logical  :: found                                    ! temporary for error check
    integer  :: iindx, jindx                             ! temporary for error check
    integer  :: ier                                      ! error status
    character(len=32) :: subname = 'surfrd'              ! subroutine name
!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg)
    allocate(pctspec(begg:endg))
    allocate(pctspecB(begg:endg))

    vegxy(:,:)   = noveg
    wtxy(:,:)  = 0._r8
    pctspec(:) = 0._r8
    pctspecB(:) = 0._r8

    if (masterproc) then
       write (6,*) 'Attempting to read surface boundary data .....'
       if (lfsurdat == ' ') then
          write(6,*)'lfsurdat must be specified' 
          call endrun()
       endif
    endif

    ! Read surface data

    if (masterproc) then
       call getfil( lfsurdat, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )
    end if

    ! Obtain surface dataset special landunit info

    call surfrd_wtxy_special(ncid, domain)

    ! Obtain surface dataset vegetated landunit info

#if (! defined DGVM)     
    if (allocate_all_vegpfts) then
       call surfrd_wtxy_veg_all(ncid, domain)
    else
       call surfrd_wtxy_veg_rank(ncid, domain)
    end if
#else
    call surfrd_wtxy_veg_dgvm(domain)
#endif

    if ( masterproc )then
       call check_ret(nf_close(ncid), subname)
       write (6,*) 'Successfully read surface boundary data'
       write (6,*)
    end if

  end subroutine surfrd

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_grid
!
! !INTERFACE:
  subroutine surfrd_get_grid(domain,filename)
!
! !DESCRIPTION:
! Read the surface dataset grid related information:
! o real edges of grid
! o integer  number of longitudes per latitude
! o real latitude  of grid cell (degrees)
! o real longitude of grid cell (degrees)
! For offline mode only the the grid does not have to be global.
! If grid is read in from dataset, grid is assumed to be global but
! does not have to be regular.
! If grid is generated by model, grid does not have to be global but 
! must then define the north, east, south, and west edges:
!
! o edges(1)    = northern edge of grid (degrees): >  -90 and <= 90
! o edges(2)    = eastern edge of grid (degrees) : see following notes
! o edges(3)    = southern edge of grid (degrees): >= -90 and <  90
! o edges(4)    = western edge of grid (degrees) : see following notes
!
!   For partial grids, northern and southern edges are any latitude
!   between 90 (North Pole) and -90 (South Pole). Western and eastern
!   edges are any longitude between -180 and 180, with longitudes
!   west of Greenwich negative. That is, western edge >= -180 and < 180;
!   eastern edge > western edge and <= 180.
!
!   For global grids, northern and southern edges are 90 (North Pole)
!   and -90 (South Pole). The western and eastern edges depend on
!   whether the grid starts at Dateline or Greenwich. Regardless,
!   these edges must span 360 degrees. Examples:
!
!                              West edge    East edge
!                            --------------------------------------------------
!  (1) Dateline            :        -180 to 180       (negative W of Greenwich)
!  (2) Greenwich (centered):    0 - dx/2 to 360 - dx/2
!
!    Grid 1 is the grid for offline mode
!    Grid 2 is the grid for cam and csm mode since the NCAR CAM
!    starts at Greenwich, centered on Greenwich
!
! !USES:
    use clm_varcon, only : spval
    use domainMod , only : domain_type,domain_init
    use areaMod   , only : celledge, cellarea                      
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ni,nj               ! size of grid on file
    integer :: ncid,dimid,varid    ! netCDF id's
    integer :: ier                 ! error status
    character(len=256)  :: locfn   ! local file name
    integer :: ret, time_index
    logical :: AREAset             ! true if area read from grid file
    logical :: NSEWset             ! true if lat/lon NSEW read from grid file
    logical :: EDGEset             ! true if EDGE NSEW read from grid file
    integer  :: strt3(3)           ! Start index to read in
    integer  :: cnt3(3)            ! Number of points to read in
    integer  :: strt1, cnt1        ! Start and count to read in for scalar
    real(r8) :: closelat           ! Single-column latitude value
    real(r8) :: closelon           ! Single-column longitude value
    integer  :: closelatidx        ! Single-column latitude index to retrieve
    integer  :: closelonidx        ! Single-column longitude index to retrieve
    character(len=32) :: subname = 'surfrd_get_grid'     ! subroutine name
!-----------------------------------------------------------------------

    AREAset = .false.
    NSEWset = .false.
    EDGEset = .false.

    if (masterproc) then

       if (filename == ' ') then
          write(6,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif

       call getfil( filename, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )

       if (single_column) then
          ni = lsmlon
          nj = lsmlat
       else
          call check_ret(nf_inq_dimid (ncid, 'lsmlon', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
          call check_ret(nf_inq_dimid (ncid, 'lsmlat', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)
       endif

    endif

    call mpi_bcast (ni, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (nj, 1, MPI_INTEGER, 0, mpicom, ier)

    call domain_init(domain,ni,nj)

    if (masterproc) then

       domain%edges(:) = spval
       ier = nf_inq_varid (ncid, 'EDGEN', varid)
       if (ier == NF_NOERR) then
          EDGEset = .true.
          call check_ret(nf_inq_varid(ncid, 'EDGEN', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%edges(1)), subname)
          call check_ret(nf_inq_varid(ncid, 'EDGEE', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%edges(2)), subname)
          call check_ret(nf_inq_varid(ncid, 'EDGES', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%edges(3)), subname)
          call check_ret(nf_inq_varid (ncid, 'EDGEW', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%edges(4)), subname)
          if (maxval(domain%edges) > 1.0e35) EDGEset = .false. !read garbage
       endif

       strt3(1)=1
       strt3(2)=1
       strt3(3)=1
       cnt3(1)=domain%ni
       cnt3(2)=domain%nj
       cnt3(3)=1
       strt1=1
       cnt1=domain%nj

       if (single_column) then
          call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
          strt3(1)=closelonidx
          strt3(2)=closelatidx
          strt1=closelatidx
          cnt1=1
       endif


       ier = nf_inq_varid (ncid, 'LATN', varid)
       if (ier == NF_NOERR) then
          NSEWset = .true.
          call check_ret(nf_inq_varid(ncid, 'LATN', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid,strt3,cnt3, domain%latn), subname)
          call check_ret(nf_inq_varid(ncid, 'LONE', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid,strt3,cnt3, domain%lone), subname)
          call check_ret(nf_inq_varid(ncid, 'LATS', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid,strt3,cnt3, domain%lats), subname)
          call check_ret(nf_inq_varid(ncid, 'LONW', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid,strt3,cnt3, domain%lonw), subname)
       endif

       ier = nf_inq_varid (ncid, 'AREA', varid)
       if (ier == NF_NOERR) then
          AREAset = .true.
          call check_ret(nf_inq_varid(ncid, 'AREA', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid,strt3,cnt3, domain%area), subname)
       endif

       call check_ret(nf_inq_varid(ncid, 'LONGXY' , varid), subname)
       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, domain%lonc), subname)
       
       call check_ret(nf_inq_varid(ncid, 'LATIXY', varid), subname)
       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, domain%latc), subname)
         
! set mask to 1 everywhere by default, override if LANDMASK exists
! if landmask exists, use it to set pftm (for older datasets)
! pftm should be overwritten below for newer datasets
       domain%mask = 1
       ier = nf_inq_varid(ncid, 'LANDMASK', varid)
       if (ier == NF_NOERR) then
          call check_ret(nf_get_vara_int(ncid, varid, strt3, cnt3, domain%mask), subname)
       endif

       ier = nf_inq_varid (ncid, 'PFTDATA_MASK', varid)
       if (ier == NF_NOERR) then
         call check_ret(nf_get_vara_int(ncid, varid, strt3, cnt3, domain%pftm), subname)
       endif
       domain%pftm = domain%mask
       where (domain%mask <= 0)
         domain%pftm = -1
       end where
       
!tcx fix, this or a test/abort should be added so overlaps can be computed
!tcx fix, this is demonstrated not bfb in cam bl311 test.
!tcx fix, see also lat_o_local in areaMod.F90
#if (1 == 0)
       ! Check lat limited to -90,90
       if (minval(domain%latc) < -90.0_r8 .or. &
           maxval(domain%latc) >  90.0_r8) then
           write(6,*) trim(subname),' Limiting lat/lon to [-90/90] from ', &
              minval(domain%latc),maxval(domain%latc)
           where (domain%latc < -90.0_r8) domain%latc = -90.0_r8
           where (domain%latc >  90.0_r8) domain%latc =  90.0_r8
       endif
#endif

       ! -------------------------------------------------------------------
       ! Define edges and area of grid cells
       ! -------------------------------------------------------------------

#if (defined SEQ_MCT) || (SEQ_ESMF) || (defined COUP_CSM)
       if (.not.NSEWset) call celledge (domain)
       if (.not.AREAset) call cellarea (domain)
#endif

#if (defined OFFLINE)
       if (.not.NSEWset) then    
          if (.not.EDGEset) then           ! global grid without use of edges
            call celledge (domain)
          else                             ! regional regular grid 
            call celledge (domain, &
                           domain%edges(1), domain%edges(2), &
                           domain%edges(3), domain%edges(4))
          end if	
       endif
       if (.not.AREAset) then
          call cellarea (domain)
       end if	
#endif

       call check_ret(nf_close(ncid), subname)

    end if   ! end of if-masterproc block

    call mpi_bcast (domain%latn , size(domain%latn) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%lats , size(domain%lats) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%lonw , size(domain%lonw) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%lone , size(domain%lone) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%area , size(domain%area) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%latc , size(domain%latc) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%lonc , size(domain%lonc) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%pftm , size(domain%pftm) , MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (domain%edges, size(domain%edges), MPI_REAL8  , 0, mpicom, ier)

  end subroutine surfrd_get_grid

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_latlon
!
! !INTERFACE:
  subroutine surfrd_get_latlon(latlon,filename,mask,mfilename,pftmflag)
!
! !DESCRIPTION:
! Read the surface dataset grid related information:
!
! !USES:
    use clm_varcon, only : spval
    use domainMod , only : latlon_type, latlon_init
    use areaMod   , only : celledge
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(latlon_type)         ,intent(inout) :: latlon   ! domain to init
    character(len=*)          ,intent(in)    :: filename ! grid filename
    integer,pointer  ,optional               :: mask(:)
    character(len=*) ,optional,intent(in)    :: mfilename ! grid filename
    logical          ,optional,intent(in)    :: pftmflag   ! is mask pft mask?
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ni,nj               ! size of grid on file
    integer :: n                   ! index
    integer :: ncid,dimid,varid    ! netCDF id's
    integer :: ncidm               ! mask file netCDF id's
    integer :: ier                 ! error status
    character(len=256)  :: locfn   ! local file name
    integer :: ret, time_index
    real(r8),allocatable :: rdata(:,:) ! temporary data
    logical :: NSEWset             ! true if lat/lon NSEW read from grid file
    logical :: EDGEset             ! true if EDGE NSEW read from grid file
    logical :: lpftmflag           ! is mask a pft mask, local copy
    integer  :: start2(2)          ! Start index to read in
    integer  :: count2(2)          ! Number of points to read in
    real(r8) :: closelat           ! Single-column latitude value
    real(r8) :: closelon           ! Single-column longitude value
    integer  :: closelatidx        ! Single-column latitude index to retrieve
    integer  :: closelonidx        ! Single-column longitude index to retrieve
    character(len=32) :: subname = 'surfrd_get_latlon'    ! subroutine name
!-----------------------------------------------------------------------

    NSEWset = .false.
    EDGEset = .false.
    lpftmflag = .false.
    ni = 0
    nj = 0

    if (present(pftmflag)) then
       lpftmflag = pftmflag
    endif

    if (masterproc) then

       if (filename == ' ') then
          write(6,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif

       call getfil( filename, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )

       if (single_column) then
          ni = lsmlon
          nj = lsmlat
       else
          ier = nf_inq_dimid (ncid, 'lon', dimid)
          if (ier == NF_NOERR) then
            call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
          endif
          ier = nf_inq_dimid (ncid, 'lat', dimid)
          if (ier == NF_NOERR) then
            call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)
          endif
          ier = nf_inq_dimid (ncid, 'lsmlon', dimid)
          if (ier == NF_NOERR) then
            call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
          endif
          ier = nf_inq_dimid (ncid, 'lsmlat', dimid)
          if (ier == NF_NOERR) then
            call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)
          endif
       endif

       if (ni == 0 .or. nj == 0) then
          write(6,*) trim(subname),' ERROR: ni or nj not set',ni,nj
          call endrun()
       endif

    endif

    call mpi_bcast (ni, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (nj, 1, MPI_INTEGER, 0, mpicom, ier)

    call latlon_init(latlon,ni,nj)
    if (present(mask)) then
       allocate(mask(ni*nj))
       mask = 1
    endif

    if (masterproc) then
       if(.not.allocated(rdata)) allocate(rdata(ni,nj))

       start2(1)  = 1
       start2(2)  = 1
       count2(1)  = ni
       count2(2)  = nj

       if (single_column) then
          call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
          start2(1) = closelonidx
          start2(2) = closelatidx
       endif

       call check_ret(nf_inq_varid(ncid, 'LONGXY' , varid), subname)
       call check_ret(nf_get_vara_double(ncid, varid, start2, count2, rdata), subname)
       !call check_ret(nf_get_var_double(ncid, varid, rdata), subname)
       latlon%lonc(:) = rdata(:,1)
          
       call check_ret(nf_inq_varid(ncid, 'LATIXY', varid), subname)
       call check_ret(nf_get_vara_double(ncid, varid, start2, count2, rdata), subname)
       !call check_ret(nf_get_var_double(ncid, varid, rdata), subname)
       latlon%latc(:) = rdata(1,:)

       latlon%edges(:) = spval
       ier = nf_inq_varid (ncid, 'EDGEN', varid)
       if (ier == NF_NOERR) then
          EDGEset = .true.
          call check_ret(nf_inq_varid(ncid, 'EDGEN', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, latlon%edges(1)), subname)
          call check_ret(nf_inq_varid(ncid, 'EDGEE', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, latlon%edges(2)), subname)
          call check_ret(nf_inq_varid(ncid, 'EDGES', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, latlon%edges(3)), subname)
          call check_ret(nf_inq_varid (ncid, 'EDGEW', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, latlon%edges(4)), subname)
          if (maxval(latlon%edges) > 1.0e35) EDGEset = .false. !read garbage
       endif

       ier = nf_inq_varid (ncid, 'LATN', varid)
       if (ier == NF_NOERR) then
          NSEWset = .true.
          call check_ret(nf_inq_varid(ncid, 'LATN', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid, start2, count2, rdata), subname)
          !call check_ret(nf_get_var_double(ncid, varid, rdata), subname)
          latlon%latn(:) = rdata(1,:)

          call check_ret(nf_inq_varid(ncid, 'LONE', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid, start2, count2, rdata), subname)
          !call check_ret(nf_get_var_double(ncid, varid, rdata), subname)
          latlon%lone(:) = rdata(:,1)

          call check_ret(nf_inq_varid(ncid, 'LATS', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid, start2, count2, rdata), subname)
          !call check_ret(nf_get_var_double(ncid, varid, rdata), subname)
          latlon%lats(:) = rdata(1,:)

          call check_ret(nf_inq_varid(ncid, 'LONW', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid, start2, count2, rdata), subname)
          !call check_ret(nf_get_var_double(ncid, varid, rdata), subname)
          latlon%lonw(:) = rdata(:,1)
       endif

#if (defined SEQ_MCT) || (defined SEQ_ESMF) || (defined COUP_CSM)
       if (.not.NSEWset) call celledge (latlon)
#endif

#if (defined OFFLINE)
       if (.not.NSEWset) then    
          if (.not.EDGEset) then           ! global grid without use of edges
            call celledge (latlon)
          else                             ! regional regular grid 
            call celledge (latlon, &
                           latlon%edges(1), latlon%edges(2), &
                           latlon%edges(3), latlon%edges(4))
          end if	
       endif
#endif
       if (present(mask)) then
          if (present(mfilename)) then
             if (mfilename == ' ') then
               write(6,*) trim(subname),' ERROR: mfilename must be specified '
               call endrun()
             endif

             call getfil( mfilename, locfn, 0 )
             call check_ret( nf_open(locfn, 0, ncidm), subname )
          else
             ncidm = ncid
          endif

          mask = 1
          ier = nf_inq_varid(ncidm, 'LANDMASK', varid)
          if (ier == NF_NOERR) then
             call check_ret(nf_get_vara_int(ncidm, varid, start2, count2, mask), subname)
          endif

          !--- if this is a pft mask, then modify and look for pftdata_mask array on dataset ---
          if (lpftmflag) then
             do n = 1,ni*nj
                if (mask(n) <= 0) mask(n) = -1
             enddo
             ier = nf_inq_varid (ncidm, 'PFTDATA_MASK', varid)
             if (ier == NF_NOERR) then
                call check_ret(nf_get_vara_int(ncidm, varid, start2, count2, mask), subname)
             endif
          endif

          if (present(mfilename)) then
             call check_ret(nf_close(ncidm), subname)
          endif

       endif

       deallocate(rdata)

!tcx fix, this or a test/abort should be added so overlaps can be computed
!tcx fix, this is demonstrated not bfb in cam bl311 test.
!tcx fix, see also lat_o_local in areaMod.F90
#if (1 == 0)
       ! Check lat limited to -90,90
       if (minval(latlon%latc) < -90.0_r8 .or. &
           maxval(latlon%latc) >  90.0_r8) then
           write(6,*) trim(subname),' Limiting lat/lon to [-90/90] from ', &
              minval(latlon%latc),maxval(latlon%latc)
           where (latlon%latc < -90.0_r8) latlon%latc = -90.0_r8
           where (latlon%latc >  90.0_r8) latlon%latc =  90.0_r8
       endif
#endif

       call check_ret(nf_close(ncid), subname)

    end if   ! end of if-masterproc block

    call mpi_bcast (latlon%latc , size(latlon%latc) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (latlon%lonc , size(latlon%lonc) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (latlon%lats , size(latlon%lats) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (latlon%latn , size(latlon%latn) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (latlon%lonw , size(latlon%lonw) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (latlon%lone , size(latlon%lone) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (latlon%edges, size(latlon%edges), MPI_REAL8  , 0, mpicom, ier)
    if (present(mask)) then
       call mpi_bcast(mask       , size(mask)         , MPI_INTEGER, 0, mpicom, ier)
    endif

  end subroutine surfrd_get_latlon




!!!!! abt rcm below
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rcmsurfrd
!
! !INTERFACE:
!  subroutine surfrd(lfsurdat, domain)
  subroutine rcmsurfrd(lakdat,glacdat,urbdat,soidat,pftdat,domain)
!
! !DESCRIPTION:
! Read the surface dataset and create subgrid weights.
! The model's surface dataset recognizes 5 basic land cover types within
! a grid cell: lake, wetland, urban, glacier, and vegetated. The vegetated
! portion of the grid cell is comprised of up to [maxpatch_pft] PFTs. These
! subgrid patches are read in explicitly for each grid cell. This is in
! contrast to LSMv1, where the PFTs were built implicitly from biome types.
!    o real edges of grid
!    o integer  number of longitudes per latitude
!    o real latitude  of grid cell (degrees)
!    o real longitude of grid cell (degrees)
!    o integer surface type: 0 = ocean or 1 = land
!    o integer soil color (1 to 20) for use with soil albedos
!    o real soil texture, %sand, for thermal and hydraulic properties
!    o real soil texture, %clay, for thermal and hydraulic properties
!    o real % of cell covered by lake    for use as subgrid patch
!    o real % of cell covered by wetland for use as subgrid patch
!    o real % of cell that is urban      for use as subgrid patch
!    o real % of cell that is glacier    for use as subgrid patch
!    o integer PFTs
!    o real % abundance PFTs (as a percent of vegetated area)
!
! !USES:
    use clm_varctl  , only : allocate_all_vegpfts
    use pftvarcon   , only : noveg
    use fileutils , only : getfil
    use domainMod , only : domain_type
    use mod_clm
    use mod_dynparam
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
!    integer , intent(out) :: vegxy(:,:)   ! PFT
!    real(r8), intent(out) :: wtxy(:,:)  ! subgrid weights
!    character(len=*), intent(in) :: lfsurdat               ! surf filename
    type(domain_type), optional, intent(in) :: domain ! domain associated with wtxy
    character(len=*), intent(in) :: soidat               ! surf filename
    character(len=*), intent(in) :: urbdat               ! surf filename
    character(len=*), intent(in) :: glacdat              ! surf filename
    character(len=*), intent(in) :: lakdat               ! surf filename
    character(len=*), intent(in) :: pftdat               ! surf filename
!    logical, intent(in) :: mask_only               ! if rank calc is done

!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: locfn                          ! local file name
    integer :: ncid,dimid,varid    ! netCDF id's
    integer  :: ncidurb,ncidlak,ncidglac,ncidsoi,ncidpft ! netCDF ids abt
    integer  :: begg,endg   
    logical  :: found                                    ! temporary for error check
    integer  :: iindx, jindx                             ! temporary for error check
    integer :: ier                 ! error status
    character(len=32) :: subname = 'rcmsurfrd'              ! subroutine name
!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg)
    allocate(pctspec(begg:endg))
    allocate(pctspecB(begg:endg))

    vegxy(:,:)   = noveg
    wtxy(:,:)  = 0._r8
    pctspec(:) = 0._r8
    pctspecB(:) = 0._r8

    if (masterproc) then
       write (6,*) 'Attempting to read surface boundary data .....'
       if (lakdat == ' ' .or. glacdat == ' ' .or. soidat == ' ' &
           .or. urbdat == ' ' .or. pftdat == ' ') then
          write(6,*)'surdat must be specified'; call endrun()
       endif
       endif

    ! Read surface data

    if (masterproc) then
!       call getfil( lfsurdat, locfn, 0 )
!       call check_ret( nf_open(locfn, 0, ncid), subname )

       call getfil( lakdat, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncidlak), subname )

       call getfil( soidat, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncidsoi), subname )

       call getfil( glacdat, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncidglac), subname )

       call getfil( urbdat, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncidurb), subname )

       call getfil( pftdat, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncidpft), subname )

    end if
          
    ! Obtain surface dataset special landunit info

    call rcmsurfrd_wtxy_special(ncidlak,ncidglac,ncidurb,ncidsoi, domain)
#if (defined VOC)
    call rcmsurfrd_bvocs(domain)   !abt
#endif
       
       
    ! Obtain surface dataset vegetated landunit info
#if (! defined DGVM)     
    if (allocate_all_vegpfts) then
       call surfrd_wtxy_veg_all(ncidpft, domain)
    else
       call surfrd_wtxy_veg_rank(ncidpft, domain)
    end if
#else
    call surfrd_wtxy_veg_dgvm(domain)
#endif
       


       ! *** abt added clm2bats_conversion call
       ! Subrountine used to convert clm land cover types
       ! to types compatable with RegCM
       ! Variables modified in RegCM in subroutine init_clmser/para.F
       ! see veg2d, veg2d1, satbrt, and satbrt1

       call clm2bats_conversion(ncidpft,ncidlak,ncidglac,ncidurb)  !abt



   
    if ( masterproc )then
!       call check_ret(nf_close(ncid), subname)
       call check_ret(nf_close(ncidsoi), subname)
       call check_ret(nf_close(ncidlak), subname)
       call check_ret(nf_close(ncidglac), subname)
       call check_ret(nf_close(ncidurb), subname)
       call check_ret(nf_close(ncidpft), subname)
       write (6,*) 'Successfully read surface boundary data'
       write (6,*)
    end if

    if(allocated(pctspec)) deallocate(pctspec) !abt added
    if(allocated(pctspecB)) deallocate(pctspecB) !abt added

  end subroutine rcmsurfrd


! abt rcm above

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rcmsurfrd_get_grid
!
! !INTERFACE:
  subroutine rcmsurfrd_get_grid(domain,filename)
!
! !DESCRIPTION:
! Read the surface dataset grid related information:
! o real edges of grid
! o integer  number of longitudes per latitude
! o real latitude  of grid cell (degrees)
! o real longitude of grid cell (degrees)
! For offline mode only the the grid does not have to be global.
! If grid is read in from dataset, grid is assumed to be global but
! does not have to be regular.
! If grid is generated by model, grid does not have to be global but 
! must then define the north, east, south, and west edges:
!
! o edges(1)    = northern edge of grid (degrees): >  -90 and <= 90
! o edges(2)    = eastern edge of grid (degrees) : see following notes
! o edges(3)    = southern edge of grid (degrees): >= -90 and <  90
! o edges(4)    = western edge of grid (degrees) : see following notes
!
!   For partial grids, northern and southern edges are any latitude
!   between 90 (North Pole) and -90 (South Pole). Western and eastern
!   edges are any longitude between -180 and 180, with longitudes
!   west of Greenwich negative. That is, western edge >= -180 and < 180;
!   eastern edge > western edge and <= 180.
!
!   For global grids, northern and southern edges are 90 (North Pole)
!   and -90 (South Pole). The western and eastern edges depend on
!   whether the grid starts at Dateline or Greenwich. Regardless,
!   these edges must span 360 degrees. Examples:
!
!                              West edge    East edge
!                            --------------------------------------------------
!  (1) Dateline            :        -180 to 180       (negative W of Greenwich)
!  (2) Greenwich (centered):    0 - dx/2 to 360 - dx/2
!
!    Grid 1 is the grid for offline mode
!    Grid 2 is the grid for cam and csm mode since the NCAR CAM
!    starts at Greenwich, centered on Greenwich
!
! !USES:
    use clm_varcon, only : spval
    use domainMod , only : domain_type,domain_init
    use areaMod   , only : celledge, cellarea                      
    use fileutils , only : getfil
    use clm_varsur, only : landmask,landfrac,satbrt_clm,r2cimask,init_tgb
    use clm_varsur, only : glatc,glonc
    use mod_clm
    use mod_dynparam
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ni,nj,nns,nni,nnj   ! size of grid on file
    integer :: ncid,dimid,varid    ! netCDF id's
    integer :: ier,counter         ! error status
    character(len=256)  :: locfn   ! local file name
    integer :: ret, time_index
    logical :: AREAset             ! true if area read from grid file
    logical :: NSEWset             ! true if lat/lon NSEW read from grid file
    logical :: EDGEset             ! true if EDGE NSEW read from grid file
    integer  :: strt3(3)           ! Start index to read in
    integer  :: cnt3(3)            ! Number of points to read in
    integer  :: strt1, cnt1        ! Start and count to read in for scalar
    real(r8),allocatable :: rcmfrac(:)
    logical ,allocatable :: mask_logic(:)
    real(r8) :: closelat           ! Single-column latitude value
    real(r8) :: closelon           ! Single-column longitude value
    integer  :: closelatidx        ! Single-column latitude index to retrieve
    integer  :: closelonidx        ! Single-column longitude index to retrieve
    character(len=32) :: subname = 'rcmsurfrd_get_grid'     ! subroutine name
!-----------------------------------------------------------------------

    AREAset = .false.
    NSEWset = .false.
    EDGEset = .false.

    if (masterproc) then

       if (filename == ' ') then
          write(6,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif

       call getfil( filename, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )

       if (single_column) then
          ni = lsmlon
          nj = lsmlat
       else
          call check_ret(nf_inq_dimid (ncid, 'lon', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
          call check_ret(nf_inq_dimid (ncid, 'lat', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)
       endif

       endif

    call mpi_bcast (ni, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (nj, 1, MPI_INTEGER, 0, mpicom, ier)

    call domain_init(domain,ni,nj)
    if(.not.associated(landmask)) allocate(landmask(ni,nj))  !abt
    if(.not.allocated(landfrac)) allocate(landfrac(ni,nj))  !abt


    if (masterproc) then
       if(.not.allocated(glatc)) allocate(glatc(ni*nj),glonc(ni*nj))

       domain%edges(:) = spval
       ier = nf_inq_varid (ncid, 'EDGEN', varid)
       if (ier == NF_NOERR) then
          EDGEset = .true.
          call check_ret(nf_inq_varid(ncid, 'EDGEN', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%edges(1)), subname)
          call check_ret(nf_inq_varid(ncid, 'EDGEE', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%edges(2)), subname)
          call check_ret(nf_inq_varid(ncid, 'EDGES', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%edges(3)), subname)
          call check_ret(nf_inq_varid (ncid, 'EDGEW', varid), subname)
          call check_ret(nf_get_var_double(ncid, varid, domain%edges(4)), subname)
          if (maxval(domain%edges) > 1.0e35) EDGEset = .false. !read garbage
       endif

       strt3(1)=1
       strt3(2)=1
       strt3(3)=1
       cnt3(1)=domain%ni
       cnt3(2)=domain%nj
       cnt3(3)=1
       strt1=1
       cnt1=domain%nj

       if (single_column) then
          call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
          strt3(1)=closelonidx
          strt3(2)=closelatidx
          strt1=closelatidx
          cnt1=1
       endif


       ier = nf_inq_varid (ncid, 'LATN', varid)
       if (ier == NF_NOERR) then
          NSEWset = .true.
          call check_ret(nf_inq_varid(ncid, 'LATN', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid,strt3,cnt3, domain%latn), subname)
          call check_ret(nf_inq_varid(ncid, 'LONE', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid,strt3,cnt3, domain%lone), subname)
          call check_ret(nf_inq_varid(ncid, 'LATS', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid,strt3,cnt3, domain%lats), subname)
          call check_ret(nf_inq_varid(ncid, 'LONW', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid,strt3,cnt3, domain%lonw), subname)
       endif

       ier = nf_inq_varid (ncid, 'AREA', varid)
       if (ier == NF_NOERR) then
          AREAset = .true.
          call check_ret(nf_inq_varid(ncid, 'AREA', varid), subname)
          call check_ret(nf_get_vara_double(ncid, varid,strt3,cnt3, domain%area), subname)
       endif

!rcm abt below
!       call check_ret(nf_inq_varid(ncid, 'LONGXY' , varid), subname)
!       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, domain%lonc), subname)
       
!       call check_ret(nf_inq_varid(ncid, 'LATIXY', varid), subname)
!       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, domain%latc), subname)

       do nnj = 1,nj
       do nni = 1,ni
         nns = (nnj-1)*ni + nni
         domain%lonc(nns) = r2cxlond_all(nni,nnj)
         domain%latc(nns) = r2cxlatd_all(nni,nnj)
         glonc(nns)       = r2cxlond_all(nni,nnj)
         glatc(nns)       = r2cxlatd_all(nni,nnj)
       enddo 
       enddo
         
       domain%area(:) = r2carea

       domain%edges(1) = r2cedgen
       domain%edges(2) = r2cedgee
       domain%edges(3) = r2cedges
       domain%edges(4) = r2cedgew

! set mask to 1 everywhere by default, override if LANDMASK exists
! if landmask exists, use it to set pftm (for older datasets)
! pftm should be overwritten below for newer datasets
       domain%mask = 1
!      ier = nf_inq_varid(ncid, 'LANDMASK', varid)
!       if (ier == NF_NOERR) then
!          call check_ret(nf_get_vara_int(ncid, varid, strt3, cnt3, domain%mask), subname)

       if(.not.allocated(mask_logic)) allocate(mask_logic(ni*nj))
       if(.not.allocated(rcmfrac))    allocate(rcmfrac(ni*nj))

       ier = nf_inq_varid(ncid, 'LANDFRAC', varid)
    if (ier == NF_NOERR) then
       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, rcmfrac), subname)

       counter     = 0
      
       do nnj = 1,nj
       do nni = 1,ni
          nns = (nnj-1)*ni + nni

          if(r2cimask == 2) then               !using weighted landfraction method
            landfrac(nni,nnj) = rcmfrac(nns)
            if(rcmfrac(nns) < 0.1) then    
              domain%mask(nns)  = 0.
              landmask(nni,nnj) = 0
              landfrac(nni,nnj) = 0.
            else
              domain%mask(nns)  = 1.
              landmask(nni,nnj) = 1
            endif

          elseif(r2cimask == 1) then           !using DOMAIN.INFO landmask
            if(satbrt_clm(nni,nnj) > 13.9 .and. satbrt_clm(nni,nnj) < 15.1) then 
              domain%mask(nns)  = 0.
              landmask(nni,nnj) = 0
              landfrac(nni,nnj) = 0.
            else
              domain%mask(nns)  = 1.
              landmask(nni,nnj) = 1
              landfrac(nni,nnj) = 1.
            endif

          endif
       enddo
       enddo

       domain%pftm = domain%mask
       where (domain%mask(:) <= 0)
         domain%pftm(:) = -1
       endwhere

       deallocate(rcmfrac) 

     endif  ! if ier not error

!rcm above

       ier = nf_inq_varid (ncid, 'PFTDATA_MASK', varid)
       if (ier == NF_NOERR) then
         call check_ret(nf_get_vara_int(ncid, varid, strt3, cnt3, domain%pftm), subname)
       endif
       
!tcx fix, this or a test/abort should be added so overlaps can be computed
!tcx fix, this is demonstrated not bfb in cam bl311 test.
!tcx fix, see also lat_o_local in areaMod.F90
#if (1 == 0)
       ! Check lat limited to -90,90
       if (minval(domain%latc) < -90.0_r8 .or. &
           maxval(domain%latc) >  90.0_r8) then
           write(6,*) trim(subname),' Limiting lat/lon to [-90/90] from ', &
              minval(domain%latc),maxval(domain%latc)
           where (domain%latc < -90.0_r8) domain%latc = -90.0_r8
           where (domain%latc >  90.0_r8) domain%latc =  90.0_r8
       endif
#endif


       ! -------------------------------------------------------------------
       ! Define edges and area of grid cells
       ! -------------------------------------------------------------------

#if (defined SEQ_MCT) || (SEQ_ESMF) || (defined COUP_CSM)
       if (.not.NSEWset) call celledge (domain)
       if (.not.AREAset) call cellarea (domain)
#endif

#if (defined OFFLINE)
!       if (.not.NSEWset) then    
!          if (.not.EDGEset) then           ! global grid without use of edges
!            call celledge (domain)
!          else                             ! regional regular grid 
            call celledge (domain, &
                           domain%edges(1), domain%edges(2), &
                           domain%edges(3), domain%edges(4))
!          end if
!       endif
!       if (.not.AREAset) then
!          call cellarea (domain)
!       end if
#endif

       call check_ret(nf_close(ncid), subname)

    end if   ! end of if-masterproc block

!abt below
    call mpi_bcast (landmask    , size(landmask)    , MPI_INTEGER, 0, mpicom, ier)  
    call mpi_bcast (landfrac    , size(landfrac)    , MPI_REAL8  , 0, mpicom, ier)  
    call mpi_bcast (satbrt_clm  , size(satbrt_clm)  , MPI_REAL8  , 0, mpicom, ier)  
    call mpi_bcast (init_tgb    , size(init_tgb)    , MPI_REAL8  , 0, mpicom, ier)  
!abt above
    call mpi_bcast (domain%latn , size(domain%latn) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%lats , size(domain%lats) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%lonw , size(domain%lonw) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%lone , size(domain%lone) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%area , size(domain%area) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%latc , size(domain%latc) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%lonc , size(domain%lonc) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (domain%pftm , size(domain%pftm) , MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (domain%edges, size(domain%edges), MPI_REAL8  , 0, mpicom, ier)


  end subroutine rcmsurfrd_get_grid





!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rcmsurfrd_get_latlon
!
! !INTERFACE:
  subroutine rcmsurfrd_get_latlon(latlon,filename,mask,mfilename,pftmflag)
!
! !DESCRIPTION:
! Read the surface dataset grid related information:
!
! !USES:
    use clm_varcon, only : spval
    use domainMod , only : latlon_type, latlon_init
    use areaMod   , only : celledge
    use fileutils , only : getfil
!abt below
    use clm_varsur, only : satbrt_clm,r2cimask
!abt above
    use mod_clm
    use mod_dynparam
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(latlon_type)         ,intent(inout) :: latlon   ! domain to init
    character(len=*)          ,intent(in)    :: filename ! grid filename
    integer,pointer  ,optional               :: mask(:)
    character(len=*) ,optional,intent(in)    :: mfilename ! grid filename
    logical          ,optional,intent(in)    :: pftmflag   ! is mask pft mask?
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: ni,nj               ! size of grid on file
    integer :: n,nni,nnj,nns        ! index
    integer :: ncid,dimid,varid    ! netCDF id's
    integer :: ncidm               ! mask file netCDF id's
    integer :: ier                 ! error status
    character(len=256)  :: locfn   ! local file name
    integer :: ret, time_index
    integer  :: strt3(3)           ! Start index to read in
    integer  :: cnt3(3)            ! Number of points to read in
    real(r8),allocatable :: rdata(:) ! temporary data
    real(r8),allocatable :: rdata1(:) ! temporary data
    integer :: counter
    logical :: NSEWset             ! true if lat/lon NSEW read from grid file
    logical :: EDGEset             ! true if EDGE NSEW read from grid file
    logical :: lpftmflag           ! is mask a pft mask, local copy
    integer  :: start2(2)          ! Start index to read in
    real(r8),allocatable :: rcmfrac(:)         ! landfraction read in from RCMnavy abt
    logical ,allocatable :: mask_logic(:)
    integer  :: count2(2)          ! Number of points to read in
    real(r8) :: closelat           ! Single-column latitude value
    real(r8) :: closelon           ! Single-column longitude value
    integer  :: closelatidx        ! Single-column latitude index to retrieve
    integer  :: closelonidx        ! Single-column longitude index to retrieve
    character(len=32) :: subname = 'rcmsurfrd_get_latlon'    ! subroutine name
!-----------------------------------------------------------------------

    NSEWset = .false.
    EDGEset = .false.
    lpftmflag = .false.
    ni = 0
    nj = 0

    if (present(pftmflag)) then
       lpftmflag = pftmflag
    endif

    if (masterproc) then

       if (filename == ' ') then
          write(6,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif

       call getfil( filename, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )

       if (single_column) then
          ni = lsmlon
          nj = lsmlat
       else
          ier = nf_inq_dimid (ncid, 'lon', dimid)
          if (ier == NF_NOERR) then
            call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
          endif
          ier = nf_inq_dimid (ncid, 'lat', dimid)
          if (ier == NF_NOERR) then
            call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)
          endif
          ier = nf_inq_dimid (ncid, 'lsmlon', dimid)
          if (ier == NF_NOERR) then
            call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
          endif
          ier = nf_inq_dimid (ncid, 'lsmlat', dimid)
          if (ier == NF_NOERR) then
            call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)
          endif
       endif

       if (ni == 0 .or. nj == 0) then
          write(6,*) trim(subname),' ERROR: ni or nj not set',ni,nj
          call endrun()
       endif

    endif

    call mpi_bcast (ni, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (nj, 1, MPI_INTEGER, 0, mpicom, ier)

    call latlon_init(latlon,ni,nj)
    if (present(mask)) then
       allocate(mask(ni*nj))
       mask = 1
    endif

    if (masterproc) then
       if(.not.allocated(rdata)) allocate(rdata(ni))
       if(.not.allocated(rdata1)) allocate(rdata1(nj))


       start2(1)  = 1
       start2(2)  = 1
       count2(1)  = ni
       count2(2)  = nj

       if (single_column) then
          call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
          start2(1) = closelonidx
          start2(2) = closelatidx
       endif

!abt rcm below       call check_ret(nf_inq_varid(ncid, 'lon', varid), subname)
!       call check_ret(nf_get_vara_double(ncid, varid, start2, count2, rdata), subname)
!       !call check_ret(nf_get_var_double(ncid, varid, rdata), subname)
!      latlon%lonc(:) = rdata(:)


!       call check_ret(nf_inq_varid(ncid, 'lat', varid), subname)
!       call check_ret(nf_get_vara_double(ncid, varid, rdata1), subname)
!       call check_ret(nf_get_var_double(ncid, varid, rdata1), subname)
!     latlon%latc(:) = rdata1(:)


       latlon%latc(:) = r2cxlatd_all(1,1:nj)
       latlon%lonc(:) = r2cxlond_all(1:ni,1)


!       latlon%edges(:) = spval
!       ier = nf_inq_varid (ncid, 'EDGEN', varid)
!       if (ier == NF_NOERR) then
!          EDGEset = .true.
!          call check_ret(nf_inq_varid(ncid, 'EDGEN', varid), subname)
!          call check_ret(nf_get_var_double(ncid, varid, latlon%edges(1)), subname)
!          call check_ret(nf_inq_varid(ncid, 'EDGEE', varid), subname)
!          call check_ret(nf_get_var_double(ncid, varid, latlon%edges(2)), subname)
!          call check_ret(nf_inq_varid(ncid, 'EDGES', varid), subname)
!          call check_ret(nf_get_var_double(ncid, varid, latlon%edges(3)), subname)
!          call check_ret(nf_inq_varid (ncid, 'EDGEW', varid), subname)
!          call check_ret(nf_get_var_double(ncid, varid, latlon%edges(4)), subname)
!          if (maxval(latlon%edges) > 1.0e35) EDGEset = .false. !read garbage
!       endif
! abt above
! abt rcm edges used 
       latlon%edges(:) = spval
       latlon%edges(1) = r2cedgen
       latlon%edges(2) = r2cedgee
       latlon%edges(3) = r2cedges
       latlon%edges(4) = r2cedgew

! abt below      ier = nf_inq_varid (ncid, 'LATN', varid)
!       if (ier == NF_NOERR) then
!          NSEWset = .true.
!          call check_ret(nf_inq_varid(ncid, 'LATN', varid), subname)
!          call check_ret(nf_get_vara_double(ncid, varid, start2, count2, rdata), subname)
!          !call check_ret(nf_get_var_double(ncid, varid, rdata), subname)
!          latlon%latn(:) = rdata(1,:)
!
!          call check_ret(nf_inq_varid(ncid, 'LONE', varid), subname)
!          call check_ret(nf_get_vara_double(ncid, varid, start2, count2, rdata), subname)
!          !call check_ret(nf_get_var_double(ncid, varid, rdata), subname)
!          latlon%lone(:) = rdata(:,1)

!          call check_ret(nf_inq_varid(ncid, 'LATS', varid), subname)
!          call check_ret(nf_get_vara_double(ncid, varid, start2, count2, rdata), subname)
!          !call check_ret(nf_get_var_double(ncid, varid, rdata), subname)
!          latlon%lats(:) = rdata(1,:)

!          call check_ret(nf_inq_varid(ncid, 'LONW', varid), subname)
!          call check_ret(nf_get_vara_double(ncid, varid, start2, count2, rdata), subname)
!          !call check_ret(nf_get_var_double(ncid, varid, rdata), subname)
!          latlon%lonw(:) = rdata(:,1)
! abt above       endif

#if (defined SEQ_MCT) || (defined SEQ_ESMF) || (defined COUP_CSM)
       if (.not.NSEWset) call celledge (latlon)
#endif

#if (defined OFFLINE)
!abt       if (.not.NSEWset) then    
!          if (.not.EDGEset) then           ! global grid without use of edges
!            call celledge (latlon)
!abt          else                             ! regional regular grid 
            call celledge (latlon, &
                           latlon%edges(1), latlon%edges(2), &
                           latlon%edges(3), latlon%edges(4))
!abt          end if
!abt       endif
#endif
       if (present(mask)) then
          if (present(mfilename)) then
             if (mfilename == ' ') then
               write(6,*) trim(subname),' ERROR: mfilename must be specified '
               call endrun()
             endif

             call getfil( mfilename, locfn, 0 )
             call check_ret( nf_open(locfn, 0, ncidm), subname )
          else
             ncidm = ncid
          endif

!abt          ier = nf_inq_varid(ncidm, 'LANDMASK', varid)
!          if (ier == NF_NOERR) then
!             call check_ret(nf_get_vara_int(ncidm, varid, start2, count2, mask), subname)
!abt          endif

           strt3(1)=1
           strt3(2)=1
           strt3(3)=1
           cnt3(1)=ni
           cnt3(2)=nj
           cnt3(3)=1

           if(.not.allocated(rcmfrac))    allocate(rcmfrac(ni*nj))
           if(.not.allocated(mask_logic)) allocate(mask_logic(ni*nj))
           ier = nf_inq_varid(ncid, 'LANDFRAC', varid)
         if (ier == NF_NOERR) then
            call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, rcmfrac), subname)


           mask = 1
           do nnj = 1,nj
           do nni = 1,ni
              nns = (nnj-1)*ni + nni

              if(r2cimask == 2) then               !using weighted landfraction method
                if(rcmfrac(nns) < 0.1) then
                  mask(nns)  = 0.
                else
                  mask(nns)  = 1.
                endif

              elseif(r2cimask == 1) then           !using DOMAIN.INFO landmask
                if(satbrt_clm(nni,nnj) > 13.9 .and. satbrt_clm(nni,nnj) < 15.1) then
                  mask(nns)  = 0.
                else
                  mask(nns)  = 1.
                endif

              endif

           enddo
           enddo


       
         endif  ! if ier not error
!rcm above


          !--- if this is a pft mask, then modify and look for pftdata_mask array on dataset ---
          if (lpftmflag) then
             do n = 1,ni*nj
                if (mask(n) <= 0) mask(n) = -1
             enddo
             ier = nf_inq_varid (ncidm, 'PFTDATA_MASK', varid)
             if (ier == NF_NOERR) then
                call check_ret(nf_get_vara_int(ncidm, varid, start2, count2, mask), subname)
             endif
          endif

          if (present(mfilename)) then
             call check_ret(nf_close(ncidm), subname)
          endif

         deallocate(rcmfrac)

       endif   !if mask present

       deallocate(rdata)
       deallocate(rdata1)

!tcx fix, this or a test/abort should be added so overlaps can be computed
!tcx fix, this is demonstrated not bfb in cam bl311 test.
!tcx fix, see also lat_o_local in areaMod.F90
#if (1 == 0)
       ! Check lat limited to -90,90
       if (minval(latlon%latc) < -90.0_r8 .or. &
           maxval(latlon%latc) >  90.0_r8) then
           write(6,*) trim(subname),' Limiting lat/lon to [-90/90] from ', &
              minval(latlon%latc),maxval(latlon%latc)
           where (latlon%latc < -90.0_r8) latlon%latc = -90.0_r8
           where (latlon%latc >  90.0_r8) latlon%latc =  90.0_r8
       endif
#endif

       call check_ret(nf_close(ncid), subname)


    end if   ! end of if-masterproc block


    call mpi_bcast (latlon%latc , size(latlon%latc) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (latlon%lonc , size(latlon%lonc) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (latlon%lats , size(latlon%lats) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (latlon%latn , size(latlon%latn) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (latlon%lonw , size(latlon%lonw) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (latlon%lone , size(latlon%lone) , MPI_REAL8  , 0, mpicom, ier)
    call mpi_bcast (latlon%edges, size(latlon%edges), MPI_REAL8  , 0, mpicom, ier)
    if (present(mask)) then
       call mpi_bcast(mask       , size(mask)         , MPI_INTEGER, 0, mpicom, ier)
    endif

  end subroutine rcmsurfrd_get_latlon


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rcmsurfrd_get_frac
!
! !INTERFACE:
  subroutine rcmsurfrd_get_frac(domain,filename)
!
! !DESCRIPTION:
! Read the landfrac dataset grid related information:
! Assume domain has already been initialized and read
!
! !USES:
    use domainMod , only : domain_type
    use fileutils , only : getfil
!abt below
    use clm_varsur, only : satbrt_clm,r2cimask
!abt above
    use mod_clm
    use mod_dynparam
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n,nni,nnj,nns       ! indices
    integer :: ni,nj,ns            ! size of grid on file
    integer :: ncid,dimid,varid    ! netCDF id's
    integer :: ier                 ! error status
    real(r8):: eps = 1.0e-12_r8    ! lat/lon error tolerance
    real(r8),allocatable:: lonc(:),latc(:)  ! local lat/lon
    character(len=256)  :: locfn   ! local file name
    integer :: ret, time_index
    integer  :: strt3(3)           ! Start index to read in
    integer  :: cnt3(3)            ! Number of points to read in
    integer  :: strt1, cnt1        ! Start and count to read in for scalar
    integer  :: counter
    logical,allocatable :: mask_logic(:)
    real(r8) :: closelat           ! Single-column latitude value
    real(r8) :: closelon           ! Single-column longitude value
    integer  :: closelatidx        ! Single-column latitude index to retrieve
    integer  :: closelonidx        ! Single-column longitude index to retrieve
    character(len=32) :: subname = 'rcmsurfrd_get_frac'     ! subroutine name
!-----------------------------------------------------------------------
 

    if (masterproc) then

       if (filename == ' ') then
          write(6,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif

       call getfil( filename, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )

       if (single_column) then
          ni = lsmlon
          nj = lsmlat
       else
          call check_ret(nf_inq_dimid (ncid, 'lon', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
          call check_ret(nf_inq_dimid (ncid, 'lat', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)
       endif

       ns = ni*nj

       if (domain%ni /= ni .or. domain%nj /= nj .or. domain%ns /= ns) then
          write(6,*) trim(subname),' ERROR: landfrac file mismatch ni,nj',domain%ni,ni,domain%nj,nj,domain%ns,ns
          call endrun()
       endif

       if(.not.allocated(latc)) allocate(latc(ni*nj),lonc(ni*nj))

       strt3(1)=1
       strt3(2)=1
       strt3(3)=1
       cnt3(1)=domain%ni
       cnt3(2)=domain%nj
       cnt3(3)=1
       strt1=1
       cnt1=domain%nj

       if (single_column) then
          call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
          strt3(1)=closelonidx
          strt3(2)=closelatidx
          strt1=closelatidx
          cnt1=1
       endif

!abt       call check_ret(nf_inq_varid(ncid, 'LONGXY' , varid), subname)
!abt       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, lonc), subname)
     
!abt     call check_ret(nf_inq_varid(ncid, 'LATIXY', varid), subname)
!abt       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, latc), subname)

!rcm abt below
       do nnj = 1,nj
       do nni = 1,ni
         nns = (nnj-1)*ni + nni
         lonc(nns) = r2cxlond_all(nni,nnj)
         latc(nns) = r2cxlatd_all(nni,nnj)
       enddo 
       enddo
!rcm above         

       do n = 1,ns
          if (abs(latc(n)-domain%latc(n)) > eps .or. &
               abs(lonc(n)-domain%lonc(n)) > eps) then
             write(6,*) trim(subname),' ERROR: landfrac file mismatch lat,lon',latc(n),domain%latc(n),lonc(n),domain%lonc(n),eps
             call endrun()
          endif
       enddo
       
!abt below       call check_ret(nf_inq_varid(ncid, 'LANDMASK', varid), subname)
!       call check_ret(nf_get_vara_int(ncid, varid, strt3, cnt3, domain%mask), subname)
       
       call check_ret(nf_inq_varid(ncid, 'LANDFRAC', varid), subname)
       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, domain%frac), subname)
 
       do nnj = 1,nj
       do nni = 1,ni
          nns = (nnj-1)*ni + nni

          if(r2cimask == 2) then               !using weighted landfraction method
            if(domain%frac(nns) < 0.1) then
              domain%mask(nns)  = 0.
            else
              domain%mask(nns)  = 1.
            endif

          elseif(r2cimask == 1) then           !using DOMAIN.INFO landmask
            if(satbrt_clm(nni,nnj) > 13.9 .and. satbrt_clm(nni,nnj) < 15.1) then
              domain%mask(nns)  = 0.
            else
              domain%mask(nns)  = 1.
            endif

          endif

       enddo
       enddo

       if(r2cimask == 2) then
         where(domain%mask(:) == 0.)
           domain%frac(:) = 0.
         endwhere
       elseif(r2cimask == 1) then           !using DOMAIN.INFO landmask
         where(domain%mask(:) == 0.)
           domain%frac(:) = 0.
         else where
           domain%frac(:) = 1.
         endwhere
       end if
!abt above

       deallocate(latc,lonc)

       call check_ret(nf_close(ncid), subname)

    end if   ! end of if-masterproc block

    call mpi_bcast (domain%mask , size(domain%mask) , MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (domain%frac , size(domain%frac) , MPI_REAL8  , 0, mpicom, ier)

  end subroutine rcmsurfrd_get_frac


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rcmsurfrd_get_topo
!
! !INTERFACE:
  subroutine rcmsurfrd_get_topo(domain,filename)
!
! !DESCRIPTION:
! Read the topo dataset grid related information:
! Assume domain has already been initialized and read
!
! !USES:
    use domainMod , only : domain_type
    use fileutils , only : getfil
    use clm_varcon, only : grav
    use clm_varsur, only : ht_rcm
    use mod_clm
    use mod_dynparam
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n,nni,nnj,nns                   ! indices abt
    integer :: ni,nj,ns            ! size of grid on file
    integer :: ncid,dimid,varid    ! netCDF id's
    integer :: ier                 ! error status
    real(r8):: eps = 1.0e-12_r8    ! lat/lon error tolerance
    real(r8),allocatable:: lonc(:),latc(:)      ! local lat/lon
    real(r8),allocatable:: TOPO(:,:)        !abt
    character(len=256)  :: locfn   ! local file name
    integer :: ret, time_index
    integer  :: strt3(3)           ! Start index to read in
    integer  :: cnt3(3)            ! Number of points to read in
    integer  :: strt1, cnt1        ! Start and count to read in for scalar
    real(r8) :: closelat           ! Single-column latitude value
    real(r8) :: closelon           ! Single-column longitude value
    integer  :: closelatidx        ! Single-column latitude index to retrieve
    integer  :: closelonidx        ! Single-column longitude index to retrieve
    character(len=32) :: subname = 'rcmsurfrd_get_topo'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then

       if (filename == ' ') then
          write(6,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif

       call getfil( filename, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )

       if (single_column) then
          ni = lsmlon
          nj = lsmlat
       else
          call check_ret(nf_inq_dimid (ncid, 'lon', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
          call check_ret(nf_inq_dimid (ncid, 'lat', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)
       endif

       ns = ni*nj

       if (domain%ni /= ni .or. domain%nj /= nj .or. domain%ns /= ns) then
          write(6,*) trim(subname),' ERROR: topo file mismatch ni,nj',domain%ni,ni,domain%nj,nj,domain%ns,ns
          call endrun()
       endif

       if(.not.allocated(latc)) allocate(latc(ns),lonc(ns))

       strt3(1)=1
       strt3(2)=1
       strt3(3)=1
       cnt3(1)=domain%ni
       cnt3(2)=domain%nj
       cnt3(3)=1
       strt1=1
       cnt1=domain%nj

       if (single_column) then
          call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
          strt3(1)=closelonidx
          strt3(2)=closelatidx
          strt1=closelatidx
          cnt1=1
       endif

!abt       call check_ret(nf_inq_varid(ncid, 'LONGXY' , varid), subname)
!abt       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, lonc), subname)

!abt       call check_ret(nf_inq_varid(ncid, 'LATIXY', varid), subname)
!abt       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, latc), subname)

!rcm abt below
       do nnj = 1,nj
       do nni = 1,ni
         nns = (nnj-1)*ni + nni
         lonc(nns) = r2cxlond_all(nni,nnj)
         latc(nns) = r2cxlatd_all(nni,nnj)
       enddo 
       enddo
!rcm above


       do n = 1,ns
          if (abs(latc(n)-domain%latc(n)) > eps .or. &
              abs(lonc(n)-domain%lonc(n)) > eps) then
             write(6,*) trim(subname),' ERROR: topo file mismatch lat,lon',latc(n),domain%latc(n),lonc(n),domain%lonc(n),eps
             call endrun()
          endif
       enddo

!abt       call check_ret(nf_inq_varid(ncid, 'TOPO', varid), subname)
!abt       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, domain%topo), subname)

!abt added REGCM elevation for continuity between regcm and clm
       do nnj = 1,nj
       do nni = 1,ni
         nns = (nnj-1)*ni + nni
         domain%topo(nns) = ht_rcm(nni,nnj)/grav     !convert from geopotential ht to ht
       enddo
       enddo

!       deallocate(ht_rcm)   !abt
       deallocate(latc,lonc)

    end if   ! end of if-masterproc block

    call mpi_bcast (domain%topo , size(domain%topo) , MPI_REAL8  , 0, mpicom, ier)

  end subroutine rcmsurfrd_get_topo



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rcmsurfrd_wtxy_special 
!
! !INTERFACE:
!  subroutine surfrd_wtxy_special(ncid, pctspec, vegxy, wtxy, domain)
  subroutine rcmsurfrd_wtxy_special(ncidlak,ncidglac,ncidurb,ncidsoi, domain)
!
! !DESCRIPTION:
! Determine weight with respect to gridcell of all special "pfts" as well
! as soil color and percent sand and clay
!
! !USES:
    use pftvarcon   , only : noveg
    use domainMod   , only : domain_type,ldomain
!abt below
    use clm_varsur  , only : landfrac
    use decompMod   , only : adecomp
!abt above
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
!    integer , intent(in)    :: ncid      ! netcdf file id 
    integer , intent(in)    :: ncidlak      ! netcdf file id for lake/wetland
    integer , intent(in)    :: ncidglac     ! netcdf file id for glacier
    integer , intent(in)    :: ncidurb      ! netcdf file id for urban
    integer , intent(in)    :: ncidsoi      ! netcdf file id for soil levels
!    real(r8), intent(inout) :: pctspec(:)! percent wrt gcell special lunits
!    integer , intent(inout) :: vegxy(:,:)  ! PFT
!    real(r8), intent(inout) :: wtxy(:,:) ! subgrid weights
    type(domain_type),intent(in) :: domain ! domain associated with wtxy
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: n,nl,ns                    ! indices
    integer  :: nni,nnj                    ! indices
    integer  :: begg,endg                  ! gcell beg/end
    integer  :: dimid,varid                ! netCDF id's
    integer  :: ret, time_index
    real(r8) :: nlevsoidata(nlevsoi)
    logical  :: found                      ! temporary for error check
    integer  :: nindx                      ! temporary for error check
    integer  :: ier                        ! error status
    real(r8),pointer :: pctgla(:)      ! percent of grid cell is glacier
    real(r8),pointer :: pctlak(:)      ! percent of grid cell is lake
    real(r8),pointer :: pctwet(:)      ! percent of grid cell is wetland
    real(r8),pointer :: pcturb(:)      ! percent of grid cell is urbanized
    integer  :: strt3(3)               ! Start index to read in
    integer  :: cnt3(3)                ! Number of points to read in
    real(r8) :: closelat               ! Single-column latitude value
    real(r8) :: closelon               ! Single-column longitude value
    integer  :: closelatidx            ! Single-column latitude index to retrieve
    integer  :: closelonidx            ! Single-column longitude index to retrieve
    character(len=32) :: subname = 'rcmsurfrd_wtxy_special'  ! subroutine name
!!-----------------------------------------------------------------------

    ns = domain%ns
    call get_proc_bounds(begg,endg)

    allocate(pctgla(begg:endg),pctlak(begg:endg))
    allocate(pctwet(begg:endg),pcturb(begg:endg))

    if (masterproc) then
!       call check_dim(ncid, 'nlevsoi', nlevsoi)
       call check_dim(ncidsoi, 'level', nlevsoi) ! abt
 
    end if   ! end of if-masterproc

    ! Obtain non-grid surface properties of surface dataset other than percent pft

    strt3(1)=1
    strt3(2)=1
    strt3(3)=1
    cnt3(1)=domain%ni
    cnt3(2)=domain%nj
    cnt3(3)=1


!abt    if (single_column) then
!abt       call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
!abt       strt3(1)=closelonidx
!abt       strt3(2)=closelatidx
!abt    endif
! abt rcm below
    call ncd_iolocal(ncidlak, 'PCT_WETLAND', 'read', pctwet, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    call ncd_iolocal(ncidlak, 'PCT_LAKE'   , 'read', pctlak, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    call ncd_iolocal(ncidglac, 'PCT_GLACIER', 'read', pctgla, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    call ncd_iolocal(ncidurb, 'PCT_URBAN'  , 'read', pcturb, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )


!abt set all special unit values to zero if less than 5% cover
!

    do nl = begg,endg
      nni = adecomp%gdc2i(nl) 
      nnj = adecomp%gdc2j(nl) 

      if(ldomain%mask(nl) == 0) then
         pctlak(nl) = 0._r8
         pctwet(nl) = 0._r8
         pctgla(nl) = 0._r8
         pcturb(nl) = 0._r8
      endif

      if(landfrac(nni,nnj) .eq. 0) then
         write(*,*)"landfrac =  ",landfrac(nni,nnj)," at ",nl 
         call endrun()
      elseif((landfrac(nni,nnj)*100._r8).lt.pctspecB(nl)) then
         pcturb(nl) = 0._r8
         pctgla(nl) = 0._r8
         pctlak(nl) = 0._r8
         pctwet(nl) = 100._r8   !set to 100% for the land portion of the grid
         write(*,*)"****** wetland is 100% *********"  !buggin
      else
         pctlak(nl) = pctlak(nl)/(landfrac(nni,nnj)*100._r8)
         pctwet(nl) = pctwet(nl)/(landfrac(nni,nnj)*100._r8)
         pctgla(nl) = pctgla(nl)/(landfrac(nni,nnj)*100._r8)
         pcturb(nl) = pcturb(nl)/(landfrac(nni,nnj)*100._r8)
      endif       

      if (pctlak(nl) < 5._r8) pctlak(nl) = 0._r8
      if (pctwet(nl) < 5._r8) pctwet(nl) = 0._r8
      if (pctgla(nl) < 5._r8) pctgla(nl) = 0._r8
      if (pcturb(nl) < 5._r8) pcturb(nl) = 0._r8

      pctspecB(nl) = pctwet(nl) + pctlak(nl) + pctgla(nl) + pcturb(nl)
    
    enddo

!abt rcm above

    pctspec = pctwet + pctlak + pctgla + pcturb

    ! Error check: glacier, lake, wetland, urban sum must be less than 100

    found = .false.
    do nl = begg,endg
       if (pctspec(nl) > 100._r8+1.e-04_r8) then
          found = .true.
          nindx = nl
          exit
       end if
       if (found) exit
    end do
    if ( found ) then
       write(6,*)'surfrd error: PFT cover>100 for nl=',nindx
       write(6,*)'pctspec/pctwet/pctlak/pcturb/pctgla  = ',pctspec(nl),pctwet(nl),pctlak(nl),pcturb(nl),pctgla(nl)
       write(6,*)'begg/endg = ',begg,endg
       call endrun()
    end if

    ! Error check that urban parameterization is not yet finished

    found = .false.
    do nl = begg,endg
       if (pcturb(nl) /= 0._r8) then
          found = .true.
          nindx = nl
          exit
       end if
       if (found) exit
    end do
    if ( found ) then
       write (6,*)'surfrd error: urban parameterization not implemented at nl= ',nindx,pcturb(nl)
       call endrun()
    end if

    ! Determine veg and wtxy for special landunits

    do nl = begg,endg
       vegxy(nl,npatch_urban)  = noveg
       wtxy(nl,npatch_urban)   = pcturb(nl)/100._r8

       vegxy(nl,npatch_lake)   = noveg
       wtxy(nl,npatch_lake)    = pctlak(nl)/100._r8

       vegxy(nl,npatch_wet)    = noveg
       wtxy(nl,npatch_wet)     = pctwet(nl)/100._r8

       vegxy(nl,npatch_glacier)= noveg
       wtxy(nl,npatch_glacier) = pctgla(nl)/100._r8
    end do

   deallocate(pctgla,pctlak,pctwet,pcturb)

  end subroutine rcmsurfrd_wtxy_special



! abt rcm above

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_frac
!
! !INTERFACE:
  subroutine surfrd_get_frac(domain,filename)
!
! !DESCRIPTION:
! Read the landfrac dataset grid related information:
! Assume domain has already been initialized and read
!
! !USES:
    use domainMod , only : domain_type
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n                   ! indices
    integer :: ni,nj,ns            ! size of grid on file
    integer :: ncid,dimid,varid    ! netCDF id's
    integer :: ier                 ! error status
    real(r8):: eps = 1.0e-12_r8    ! lat/lon error tolerance
    real(r8),allocatable:: lonc(:),latc(:)  ! local lat/lon
    character(len=256)  :: locfn   ! local file name
    integer :: ret, time_index
    integer  :: strt3(3)           ! Start index to read in
    integer  :: cnt3(3)            ! Number of points to read in
    integer  :: strt1, cnt1        ! Start and count to read in for scalar
    real(r8) :: closelat           ! Single-column latitude value
    real(r8) :: closelon           ! Single-column longitude value
    integer  :: closelatidx        ! Single-column latitude index to retrieve
    integer  :: closelonidx        ! Single-column longitude index to retrieve
    character(len=32) :: subname = 'surfrd_get_frac'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then

       if (filename == ' ') then
          write(6,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif

       call getfil( filename, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )

       if (single_column) then
          ni = lsmlon
          nj = lsmlat
       else
          call check_ret(nf_inq_dimid (ncid, 'lsmlon', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
          call check_ret(nf_inq_dimid (ncid, 'lsmlat', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)
       endif

       ns = ni*nj

       if (domain%ni /= ni .or. domain%nj /= nj .or. domain%ns /= ns) then
          write(6,*) trim(subname),' ERROR: landfrac file mismatch ni,nj',domain%ni,ni,domain%nj,nj,domain%ns,ns
          call endrun()
       endif

       if(.not.allocated(latc)) allocate(latc(ni*nj),lonc(ni*nj))

       strt3(1)=1
       strt3(2)=1
       strt3(3)=1
       cnt3(1)=domain%ni
       cnt3(2)=domain%nj
       cnt3(3)=1
       strt1=1
       cnt1=domain%nj

       if (single_column) then
          call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
          strt3(1)=closelonidx
          strt3(2)=closelatidx
          strt1=closelatidx
          cnt1=1
       endif

       call check_ret(nf_inq_varid(ncid, 'LONGXY' , varid), subname)
       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, lonc), subname)
          
       call check_ret(nf_inq_varid(ncid, 'LATIXY', varid), subname)
       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, latc), subname)

       do n = 1,ns
          if (abs(latc(n)-domain%latc(n)) > eps .or. &
               abs(lonc(n)-domain%lonc(n)) > eps) then
             write(6,*) trim(subname),' ERROR: landfrac file mismatch lat,lon',latc(n),domain%latc(n),lonc(n),domain%lonc(n),eps
             call endrun()
          endif
       enddo
       
       call check_ret(nf_inq_varid(ncid, 'LANDMASK', varid), subname)
       call check_ret(nf_get_vara_int(ncid, varid, strt3, cnt3, domain%mask), subname)
       
       call check_ret(nf_inq_varid(ncid, 'LANDFRAC', varid), subname)
       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, domain%frac), subname)
       

       deallocate(latc,lonc)

       call check_ret(nf_close(ncid), subname)

    end if   ! end of if-masterproc block

    call mpi_bcast (domain%mask , size(domain%mask) , MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (domain%frac , size(domain%frac) , MPI_REAL8  , 0, mpicom, ier)

  end subroutine surfrd_get_frac

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_get_topo
!
! !INTERFACE:
  subroutine surfrd_get_topo(domain,filename)
!
! !DESCRIPTION:
! Read the topo dataset grid related information:
! Assume domain has already been initialized and read
!
! !USES:
    use domainMod , only : domain_type
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    type(domain_type),intent(inout) :: domain   ! domain to init
    character(len=*) ,intent(in)    :: filename ! grid filename
!
! !CALLED FROM:
! subroutine initialize
!
! !REVISION HISTORY:
! Created by T Craig
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: n                   ! indices
    integer :: ni,nj,ns            ! size of grid on file
    integer :: ncid,dimid,varid    ! netCDF id's
    integer :: ier                 ! error status
    real(r8):: eps = 1.0e-12_r8    ! lat/lon error tolerance
    real(r8),allocatable:: lonc(:),latc(:)  ! local lat/lon
    character(len=256)  :: locfn   ! local file name
    integer :: ret, time_index
    integer  :: strt3(3)           ! Start index to read in
    integer  :: cnt3(3)            ! Number of points to read in
    integer  :: strt1, cnt1        ! Start and count to read in for scalar
    real(r8) :: closelat           ! Single-column latitude value
    real(r8) :: closelon           ! Single-column longitude value
    integer  :: closelatidx        ! Single-column latitude index to retrieve
    integer  :: closelonidx        ! Single-column longitude index to retrieve
    character(len=32) :: subname = 'surfrd_get_topo'     ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then

       if (filename == ' ') then
          write(6,*) trim(subname),' ERROR: filename must be specified '
          call endrun()
       endif

       call getfil( filename, locfn, 0 )
       call check_ret( nf_open(locfn, 0, ncid), subname )

       if (single_column) then
          ni = lsmlon
          nj = lsmlat
       else
          call check_ret(nf_inq_dimid (ncid, 'lsmlon', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, ni), subname)
          call check_ret(nf_inq_dimid (ncid, 'lsmlat', dimid), subname)
          call check_ret(nf_inq_dimlen(ncid, dimid, nj), subname)
       endif

       ns = ni*nj

       if (domain%ni /= ni .or. domain%nj /= nj .or. domain%ns /= ns) then
          write(6,*) trim(subname),' ERROR: topo file mismatch ni,nj',domain%ni,ni,domain%nj,nj,domain%ns,ns
          call endrun()
       endif

       if(.not.allocated(latc)) allocate(latc(ns),lonc(ns))

       strt3(1)=1
       strt3(2)=1
       strt3(3)=1
       cnt3(1)=domain%ni
       cnt3(2)=domain%nj
       cnt3(3)=1
       strt1=1
       cnt1=domain%nj

       if (single_column) then
          call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
          strt3(1)=closelonidx
          strt3(2)=closelatidx
          strt1=closelatidx
          cnt1=1
       endif

       call check_ret(nf_inq_varid(ncid, 'LONGXY' , varid), subname)
       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, lonc), subname)

       call check_ret(nf_inq_varid(ncid, 'LATIXY', varid), subname)
       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, latc), subname)

       do n = 1,ns
          if (abs(latc(n)-domain%latc(n)) > eps .or. &
              abs(lonc(n)-domain%lonc(n)) > eps) then
             write(6,*) trim(subname),' ERROR: topo file mismatch lat,lon',latc(n),domain%latc(n),lonc(n),domain%lonc(n),eps
             call endrun()
          endif
       enddo

       call check_ret(nf_inq_varid(ncid, 'TOPO', varid), subname)
       call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, domain%topo), subname)

       deallocate(latc,lonc)

       call check_ret(nf_close(ncid), subname)

    end if   ! end of if-masterproc block

    call mpi_bcast (domain%topo , size(domain%topo) , MPI_REAL8  , 0, mpicom, ier)

  end subroutine surfrd_get_topo

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_special 
!
! !INTERFACE:
!  subroutine surfrd_wtxy_special(ncid, pctspec, vegxy, wtxy, domain)
  subroutine surfrd_wtxy_special(ncid, domain)
!
! !DESCRIPTION:
! Determine weight with respect to gridcell of all special "pfts" as well
! as soil color and percent sand and clay
!
! !USES:
    use pftvarcon   , only : noveg
    use domainMod   , only : domain_type
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer , intent(in)    :: ncid      ! netcdf file id 
!    real(r8), intent(inout) :: pctspec(:)! percent wrt gcell special lunits
!    integer , intent(inout) :: vegxy(:,:)  ! PFT
!    real(r8), intent(inout) :: wtxy(:,:) ! subgrid weights
    type(domain_type),intent(in) :: domain ! domain associated with wtxy
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: n,nl,ns                    ! indices
    integer  :: begg,endg                  ! gcell beg/end
    integer  :: dimid,varid                ! netCDF id's
    integer  :: ret, time_index
    real(r8) :: nlevsoidata(nlevsoi)
    logical  :: found                      ! temporary for error check
    integer  :: nindx                      ! temporary for error check
    integer  :: ier                        ! error status
    real(r8),pointer :: pctgla(:)      ! percent of grid cell is glacier
    real(r8),pointer :: pctlak(:)      ! percent of grid cell is lake
    real(r8),pointer :: pctwet(:)      ! percent of grid cell is wetland
    real(r8),pointer :: pcturb(:)      ! percent of grid cell is urbanized
    integer  :: strt3(3)               ! Start index to read in
    integer  :: cnt3(3)                ! Number of points to read in
    real(r8) :: closelat               ! Single-column latitude value
    real(r8) :: closelon               ! Single-column longitude value
    integer  :: closelatidx            ! Single-column latitude index to retrieve
    integer  :: closelonidx            ! Single-column longitude index to retrieve
    character(len=32) :: subname = 'surfrd_wtxy_special'  ! subroutine name
!!-----------------------------------------------------------------------

    ns = domain%ns
    call get_proc_bounds(begg,endg)

    allocate(pctgla(begg:endg),pctlak(begg:endg))
    allocate(pctwet(begg:endg),pcturb(begg:endg))

    if (masterproc) then
       call check_dim(ncid, 'nlevsoi', nlevsoi)
    end if   ! end of if-masterproc

    ! Obtain non-grid surface properties of surface dataset other than percent pft

    strt3(1)=1
    strt3(2)=1
    strt3(3)=1
    cnt3(1)=domain%ni
    cnt3(2)=domain%nj
    cnt3(3)=1

    if (single_column) then
       call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
       strt3(1)=closelonidx
       strt3(2)=closelatidx
    endif

    call ncd_iolocal(ncid, 'PCT_WETLAND', 'read', pctwet, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    call ncd_iolocal(ncid, 'PCT_LAKE'   , 'read', pctlak, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    call ncd_iolocal(ncid, 'PCT_GLACIER', 'read', pctgla, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    call ncd_iolocal(ncid, 'PCT_URBAN'  , 'read', pcturb, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )

    pctspec = pctwet + pctlak + pctgla + pcturb

    ! Error check: glacier, lake, wetland, urban sum must be less than 100

    found = .false.
    do nl = begg,endg
       if (pctspec(nl) > 100._r8+1.e-04_r8) then
          found = .true.
          nindx = nl
          exit
       end if
       if (found) exit
    end do
    if ( found ) then
       write(6,*)'surfrd error: PFT cover>100 for nl=',nindx
       call endrun()
    end if

    ! Error check that urban parameterization is not yet finished

    found = .false.
    do nl = begg,endg
       if (pcturb(nl) /= 0._r8) then
          found = .true.
          nindx = nl
          exit
       end if
       if (found) exit
    end do
    if ( found ) then
       write (6,*)'surfrd error: urban parameterization not implemented at nl= ',nindx,pcturb(nl)
       call endrun()
    end if

    ! Determine veg and wtxy for special landunits

    do nl = begg,endg
       vegxy(nl,npatch_urban)  = noveg
       wtxy(nl,npatch_urban)   = pcturb(nl)/100._r8

       vegxy(nl,npatch_lake)   = noveg
       wtxy(nl,npatch_lake)    = pctlak(nl)/100._r8

       vegxy(nl,npatch_wet)    = noveg
       wtxy(nl,npatch_wet)     = pctwet(nl)/100._r8

       vegxy(nl,npatch_glacier)= noveg
       wtxy(nl,npatch_glacier) = pctgla(nl)/100._r8
    end do

   deallocate(pctgla,pctlak,pctwet,pcturb)

  end subroutine surfrd_wtxy_special

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rcmsurfrd_bvocs 
!
! !INTERFACE:
  subroutine rcmsurfrd_bvocs (domain)
!
! !DESCRIPTION:
! Read in biogenic emission factor maps
!
! !USES:
    use clm_varvoc
    use clm_varpar  , only : nvoc    
    use domainMod   , only : domain_type
    use clm_varctl  , only : mksrf_fisop,mksrf_fmbo,mksrf_fapin,mksrf_fbpin
    use clm_varctl  , only : mksrf_fco,mksrf_flimo,mksrf_fsabi,mksrf_fmyrc
    use clm_varctl  , only : mksrf_focim,mksrf_facar,mksrf_fomtp,mksrf_ffarn
    use clm_varctl  , only : mksrf_facto,mksrf_fmeoh,mksrf_fosqt,mksrf_fbcar
    use clm_varctl  , only : mksrf_fmeth,mksrf_fno,mksrf_facta,mksrf_fform
    use clm_varctl  , only : nsrest
    use fileutils   , only : getfil
    use clm_time_manager, only : get_step_size
    use mod_clm
    use mod_dynparam
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
!
    type(domain_type),intent(in) :: domain ! domain associated with wtxy
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Ahmed Tawfik
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: locfn        ! local file name
    integer  :: n,nl,ns,dtime          ! indices
    integer  :: begg,endg              ! gcell beg/end
    integer  :: dimid,varid            ! netCDF id's
    integer  :: ncidv(nvoc)            ! netCDF ids
    character(len=256) :: datfil(nvoc) ! filenames
    logical  :: found                  ! temporary for error check
    integer  :: nindx                  ! temporary for error check
    integer  :: ier                    ! error status
    integer  :: strt3(3)               ! Start index to read in
    integer  :: cnt3(3)                ! Number of points to read in
    character(len=32) :: subname = 'rcmsurfrd_bvocs'  ! subroutine name
    logical  :: there                  ! checking file existence
!!-----------------------------------------------------------------------

    datfil(1)   = mksrf_fisop
    datfil(2)   = mksrf_fmyrc
    datfil(3)   = mksrf_fsabi
    datfil(4)   = mksrf_flimo
    datfil(5)   = mksrf_fco
    datfil(6)   = mksrf_fmbo
    datfil(7)   = mksrf_fapin
    datfil(8)   = mksrf_fbpin
    datfil(9)   = mksrf_focim
    datfil(10)  = mksrf_facar
    datfil(11)  = mksrf_fomtp
    datfil(12)  = mksrf_ffarn
    datfil(13)  = mksrf_fbcar
    datfil(14)  = mksrf_fosqt
    datfil(15)  = mksrf_fmeoh
    datfil(16)  = mksrf_facto
    datfil(17)  = mksrf_fmeth
    datfil(18)  = mksrf_fno
    datfil(19)  = mksrf_facta
    datfil(20)  = mksrf_fform
    ns = domain%ns
    call get_proc_bounds(begg,endg)

    allocate(ef_iso(begg:endg))
    allocate(ef_mbo(begg:endg))
    allocate(ef_bpin(begg:endg))
    allocate(ef_apin(begg:endg))
    allocate(ef_myrc(begg:endg))
    allocate(ef_limo(begg:endg))
    allocate(ef_sabi(begg:endg))
    allocate(ef_acar(begg:endg))
    allocate(ef_bcar(begg:endg))
    allocate(ef_omtp(begg:endg))
    allocate(ef_farn(begg:endg))
    allocate(ef_osqt(begg:endg))
    allocate(ef_meoh(begg:endg))
    allocate(ef_meth(begg:endg))
    allocate(ef_acto(begg:endg))
    allocate(ef_no(begg:endg))
    allocate(ef_ocim(begg:endg))
    allocate(ef_co(begg:endg))
    allocate(ef_form(begg:endg))
    allocate(ef_acta(begg:endg))

    r2cefmap(:) = 1

    if (masterproc) then
      write (6,*) 'Attempting to read Biogenic Emissions data .....'
      do n = 1,nvoc
         inquire(file=trim(datfil(n)),exist=there)
         if (there) then
           call getfil( datfil(n), locfn, 0 )
           call check_ret( nf_open(locfn, 0, ncidv(n)), subname )
           write(6,*)'emission file ',trim(datfil(n)),' being used '
         else
           write(6,*)'emission file',n,' not found, using pft for emission'; r2cefmap(n)=2
         endif
      enddo 
    end if ! masterproc
    call mpi_bcast (r2cefmap    , size(r2cefmap)    , MPI_INTEGER, 0, mpicom, ier)  
    call mpi_bcast (ncidv       , size(ncidv)       , MPI_INTEGER, 0, mpicom, ier)  


    ! Obtain non-grid emission factor maps

    strt3(1)=1
    strt3(2)=1
    strt3(3)=1
    cnt3(1)=domain%ni
    cnt3(2)=domain%nj
    cnt3(3)=1

    if(r2cefmap(1) == 1) then
      call ncd_iolocal(ncidv(1), 'ISOP' , 'read', ef_iso , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(6) == 1) then
      call ncd_iolocal(ncidv(6), 'MBO'  , 'read', ef_mbo , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(8) == 1) then
      call ncd_iolocal(ncidv(8),'BPINE', 'read', ef_bpin, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(7) == 1) then
      call ncd_iolocal(ncidv(7),'APIN', 'read', ef_apin, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(2) == 1) then
      call ncd_iolocal(ncidv(2), 'MYRCENE'   , 'read', ef_myrc , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(4) == 1) then
      call ncd_iolocal(ncidv(4), 'LIMO'  , 'read', ef_limo , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(3) == 1) then
      call ncd_iolocal(ncidv(3),'SABINENE'  , 'read', ef_sabi , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(9) == 1) then
      call ncd_iolocal(ncidv(9),'OCIMENE'   , 'read', ef_ocim , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(10) == 1) then
      call ncd_iolocal(ncidv(10), 'A3CARENE'  , 'read', ef_acar , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(11) == 1) then
      call ncd_iolocal(ncidv(11), 'O_MONO'    , 'read', ef_omtp , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(12) == 1) then
      call ncd_iolocal(ncidv(12),'FARNI'     , 'read', ef_farn , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(13) == 1) then
      call ncd_iolocal(ncidv(13),'B_CARY'    , 'read', ef_bcar , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(14) == 1) then
      call ncd_iolocal(ncidv(14), 'O_SESQ'    , 'read', ef_osqt , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(15) == 1) then
      call ncd_iolocal(ncidv(15), 'MEOH'      , 'read', ef_meoh , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(16) == 1) then
      call ncd_iolocal(ncidv(16),'ACETONE'   , 'read', ef_acto , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(17) == 1) then
      call ncd_iolocal(ncidv(17),'METHANE'   , 'read', ef_meth , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(18) == 1) then
      call ncd_iolocal(ncidv(18), 'NO'        , 'read', ef_no   , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(19) == 1) then
      call ncd_iolocal(ncidv(19), 'ACET_ETH'  , 'read', ef_acta , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    end if
    if(r2cefmap(20) == 1) then
      call ncd_iolocal(ncidv(20),'FORM'      , 'read', ef_form , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif
    if(r2cefmap(5) == 1) then
      call ncd_iolocal(ncidv(5),'CO'        , 'read', ef_co   , begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, strt3, cnt3 )
    endif

    ! Initialize time accumulation variables
    dtime = r2cdtime
    n24   = 86400/dtime
    n240  = 864000/dtime
!    if(nsrest == 0) then
      c24   = 0
      c240  = 0
!    endif

    if (masterproc) then
     do n = 1,nvoc
       if(r2cefmap(n) == 1) call check_ret(nf_close(ncidv(n)), subname)
     enddo
       write (6,*) 'Successfully read biogenic emissions data'
       write (6,*)  
    end if

  end subroutine rcmsurfrd_bvocs

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_veg_rank
!
! !INTERFACE:
!  subroutine surfrd_wtxy_veg_rank(ncid, pctspec, vegxy, wtxy, domain)
  subroutine surfrd_wtxy_veg_rank(ncid, domain)
!
! !DESCRIPTION:
! Determine wtxy and veg arrays for non-dynamic landuse mode
!
! !USES:
    use clm_varctl, only : create_crop_landunit
    use pftvarcon   , only : crop, noveg
    use domainMod   , only : domain_type
    use mod_clm
    use mod_dynparam
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer , intent(in)    :: ncid        ! netcdf file id 
    type(domain_type),intent(in) :: domain ! domain associated with wtxy
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: k,m,k1,k2,n,nl,ns               ! indices
    integer  :: begg,endg                       ! beg/end gcell index
    integer  :: dimid,varid                     ! netCDF id's
    integer  :: start(4),count(4)               ! netCDF start/count arrays
    integer  :: cropcount                       ! temporary counter
    real(r8),allocatable :: sumvec(:)           ! temporary vector sum
    integer,allocatable :: nindx(:)
    logical  :: found                           ! temporary for error check
    integer  :: miss = 99999                    ! missing data indicator
    real(r8) :: wst(0:numpft)                   ! as pft at specific i, j
    integer ,allocatable :: wsti(:)             ! ranked indices largest wst values
    real(r8) :: wst_sum                         ! sum of %pft
    real(r8) :: sumpct                          ! sum of %pft over maxpatch_pft
    real(r8) :: diff                            ! the difference (wst_sum - sumpct)
    real(r8) :: rmax                            ! maximum patch cover
    integer ,allocatable :: pft(:,:)            ! PFT
    integer ,allocatable :: cft(:,:)            ! CFT
    real(r8),allocatable :: pctcft_lunit(:,:)   ! % of crop lunit area for CFTs
    real(r8),allocatable :: pctpft_lunit(:,:)   ! % of vegetated lunit area PFTs
!abt below
    real(r8),allocatable :: adj_pctpft(:,:)     ! % of crop lunit area for CFTs
    real(r8),allocatable :: adj_wst(:)          ! % of vegetated lunit area PFTs
    real(r8) :: adj_wst_sum
    integer  :: ipft
!abt above
    integer  :: ier                             ! error status
    real(r8),allocatable :: pctpft(:,:)         ! percent of vegetated gridcell area for PFTs
    real(r8),pointer :: arrayl(:)               ! local array
    integer ,pointer :: irrayg(:)               ! global array
    integer  :: ret, time_index
    real(r8),allocatable :: rmaxpatchdata(:)
    integer ,allocatable :: imaxpatchdata(:)
    real(r8) :: numpftp1data(0:numpft)         
    integer  :: strt3(3)                        ! Start index to read in
    integer  :: cnt3(3)                         ! Number of points to read in
    real(r8) :: closelat                        ! Single-column latitude value
    real(r8) :: closelon                        ! Single-column longitude value
    integer  :: closelatidx                     ! Single-column latitude index to retrieve
    integer  :: closelonidx                     ! Single-column longitude index to retrieve
    character(len=32) :: subname = 'surfrd_wtxy_veg_rank'  ! subroutine name

!-----------------------------------------------------------------------

    ns = domain%ns
    call get_proc_bounds(begg,endg)

    allocate(sumvec(begg:endg))
    allocate(nindx(begg:endg))
    allocate(cft(begg:endg,maxpatch_cft))
    allocate(pft(begg:endg,maxpatch_pft))
    allocate(pctcft_lunit(begg:endg,maxpatch_cft))
    allocate(pctpft_lunit(begg:endg,maxpatch_pft))
    allocate(pctpft(begg:endg,0:numpft))
    allocate(wsti(maxpatch_pft))

    if(.not.allocated(adj_wst)) allocate(adj_wst(0:numpft))               !abt
    if(.not.allocated(adj_pctpft)) allocate(adj_pctpft(begg:endg,0:numpft))  !abt

    if (masterproc) then
!       call check_dim(ncid, 'lsmpft', numpft+1) abt
       call check_dim(ncid, 'level', numpft+1)
    end if
    allocate(arrayl(begg:endg))
    do n = 0,numpft
       start(1) = 1
       count(1) = domain%ni
       start(2) = 1
       count(2) = domain%nj
       start(3) = n+1	 ! dataset is 1:numpft+1, not 0:numpft
       count(3) = 1
       count(4) = 1
       start(4) = 1   !abt added this to include time dimension

       if (single_column) then
          call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
          start(1)=closelonidx
          start(2)=closelatidx
       endif
       call ncd_iolocal(ncid, 'PCT_PFT', 'read', arrayl, begg, endg, gsMap_lnd_gdc2glo, perm_lnd_gdc2glo, start, count)
       pctpft(begg:endg,n) = arrayl(begg:endg)
    enddo
    deallocate(arrayl)

    ! 1. pctpft must go back to %vegetated landunit instead of %gridcell
    ! 2. pctpft bare = 100 when landmask = 1 and 100% special landunit
    ! NB: (1) and (2) do not apply to crops.
    ! For now keep all cfts (< 4 anyway) instead of 4 most dominant cfts

    do nl = begg,endg
      cft(nl,:) = 0
      pctcft_lunit(nl,:) = 0._r8
      cropcount = 0
      if (pctspecB(nl) < 99.0_r8) then
        do m = 0, numpft
          if (create_crop_landunit) then
            ! Separate crop landunit is to be created
            if (crop(m) == 1._r8 .and. pctpft(nl,m) > 0._r8) then
              cropcount = cropcount + 1
              if (cropcount > maxpatch_cft) then
                write(6,*) 'ERROR surfrdMod: cropcount>maxpatch_cft'
                call endrun()
              end if
              cft(nl,cropcount) = m
              pctcft_lunit(nl,cropcount) = &
                  pctpft(nl,m)*100._r8/(100._r8-pctspec(nl))
              pctpft(nl,m) = 0.0_r8
            else if (crop(m) == 0._r8) then
              pctpft(nl,m) = pctpft(nl,m) * 100._r8/(100._r8-pctspec(nl))
            end if
          else
            ! Separate crop landunit is not created
!abt        pctpft(nl,m) = pctpft(nl,m) * 100._r8/(100._r8-pctspec(nl))
            !input file already is %pft of vegetation (sum of pfts = 100%)
            pctpft(nl,m) = pctpft(nl,m)
          end if
        end do
      else
        pctpft(nl,1:numpft) = 0._r8
      end if
    end do

    ! Find pft and pct arrays 
    ! Save percent cover by PFT [wst] and total percent cover [wst_sum]

    do nl = begg,endg

      wst_sum = 0._r8
      sumpct = 0._r8
      do m = 0, numpft
        wst(m) = pctpft(nl,m)
        wst_sum = wst_sum + pctpft(nl,m)
      end do

!!!abt below!! Correct for gridcells that do not have any percent of pft cover
             ! Use BATS landuse types from DOMAIN.INFO and convert to the equivalent 
             ! CLM landuse/pft type.  For example, BATS Evergreen Shrub type would be
             ! Broadleaf Evergreen Shrub in CLM.

      if ( abs(sum(wst(:)) + pctspecB(nl) - 100.0_r8) > 1.0_r8 ) then 
        do m = 0,numpft
          wst(m)       = 0._r8
          pctpft(nl,m) = 0._r8
          wst_sum      = 0._r8
          adj_wst(m)   = 0._r8
          adj_pctpft(nl,m) = 0._r8
          wst_sum          = 0._r8
        enddo
        call pft_adjustment(begg,endg,ipft,nl,adj_wst, &
                            adj_wst_sum,adj_pctpft,domain)  !abt
        do m = 0,numpft
          wst(m)       = adj_wst(m)
          pctpft(nl,m) = adj_pctpft(nl,m) 
        enddo
        wst_sum      = adj_wst_sum
      endif

!!!abt rcm above

      if (domain%pftm(nl) >= 0) then

        ! Rank [wst] in ascendg order to obtain the top [maxpatch_pft] PFTs

        if ( sum(wst(:)) > 1.0_r8 ) then
          call surfrd_mkrank (numpft, wst, miss, wsti, maxpatch_pft)
        end if

        ! Fill in [pft] and [pctpft] with data for top [maxpatch_pft] PFTs.
        ! If land model grid cell is ocean, set to no PFTs.
        ! If land model grid cell is land then:
        !  1. If [pctlnd_o] = 0, there is no PFT data from the input grid.
        !     Since need land data, use bare ground.
        !  2. If [pctlnd_o] > 0, there is PFT data from the input grid but:
        !     a. use the chosen PFT so long as it is not a missing value
        !     b. missing value means no more PFTs with cover > 0
                 
        do m = 1, maxpatch_pft
          if (wsti(m) /=  miss) then
            pft(nl,m) = wsti(m)
            pctpft_lunit(nl,m) = wst(wsti(m))
          else
            pft(nl,m) = noveg
            pctpft_lunit(nl,m) = 0._r8
          end if
          sumpct = sumpct + pctpft_lunit(nl,m)
        end do
      else                               ! model grid wants ocean
        do m = 1, maxpatch_pft
          pft(nl,m) = 0
          pctpft_lunit(nl,m) = 0._r8
        end do
      end if

      ! Correct for the case of more than [maxpatch_pft] PFTs present
                
!       if (sumpct < wst_sum) then
!          diff  = wst_sum - sumpct
!          sumpct = 0._r8
!          do m = 1, maxpatch_pft
!             pctpft_lunit(nl,m) = pctpft_lunit(nl,m) + diff/maxpatch_pft
!             sumpct = sumpct + pctpft_lunit(nl,m)
!          end do
!       end if

      ! Error check: make sure have a valid PFT

      do m = 1,maxpatch_pft
        if (pft(nl,m) < 0 .or. pft(nl,m) > numpft) then
          write (6,*)'surfrd error: invalid PFT at gridcell nl=',nl,pft(nl,m)
          call endrun()
        end if
      end do

      ! As done in mksrfdatMod.F90 for other percentages, truncate pctpft to
      ! ensure that weight relative to landunit is not nonzero
      ! (i.e. a very small number such as 1e-16) where it really should be zero
      ! The following if-block is here to preserve roundoff level differences
      ! between the call to surfrd_wtxy_veg_all and surfrd_wtxy_veg_rank

      ! GRAZIANO
      ! This is done in clm2rcm now
      ! GRAZIANO

!     if (maxpatch_pft < numpft+1) then
!       do m=1,maxpatch_pft
!         pctpft_lunit(nl,m) = float(nint(pctpft_lunit(nl,m)))
!       end do
!       do m=1,maxpatch_cft
!         pctcft_lunit(nl,m) = float(nint(pctcft_lunit(nl,m)))
!       end do
!     end if
!                   
!     ! Make sure sum of PFT cover == 100 for land points. If not,
!     ! subtract excess from most dominant PFT.
!
!     rmax = -9999._r8
!     k1 = -9999
!     k2 = -9999
!     sumpct = 0._r8
!     do m = 1, maxpatch_pft
!       sumpct = sumpct + pctpft_lunit(nl,m)
!       if (pctpft_lunit(nl,m) > rmax) then
!         k1 = m
!         rmax = pctpft_lunit(nl,m)
!       end if
!     end do
!     do m = 1, maxpatch_cft
!       sumpct = sumpct + pctcft_lunit(nl,m)
!       if (pctcft_lunit(nl,m) > rmax) then
!         k2 = m
!         rmax = pctcft_lunit(nl,m)
!       end if
!     end do
!     if (k1 == -9999 .and. k2 == -9999) then
!       write(6,*)'surfrd error: largest PFT patch not found'
!       call endrun()
!     else if (domain%pftm(nl) >= 0) then
!       if (sumpct < 95 .or. sumpct > 105._r8) then
!         write(6,*)'surfrd error: sum of PFT cover =',sumpct,' at nl=',nl
!         call endrun()
!       else if (sumpct /= 100._r8 .and. k2 /= -9999) then
!         pctcft_lunit(nl,k2) = pctcft_lunit(nl,k2) - (sumpct-100._r8)
!       else if (sumpct /= 100._r8) then
!         pctpft_lunit(nl,k1) = pctpft_lunit(nl,k1) - (sumpct-100._r8)
!       end if
!     end if

      ! Error check: make sure PFTs sum to 100% cover

!     sumpct = 0._r8
!     do m = 1, maxpatch_pft
!       sumpct = sumpct + pctpft_lunit(nl,m)
!     end do
!     do m = 1, maxpatch_cft
!       sumpct = sumpct + pctcft_lunit(nl,m)
!     end do
!     if (domain%pftm(nl) >= 0) then
!       if (abs(sumpct - 100._r8) > 0.000001_r8) then
!         write(6,*)'surfrdMod error: sum(pct) over maxpatch_pft is not = 100.'
!         write(6,*)sumpct, nl
!         call endrun()
!       end if
!       if (sumpct < -0.000001_r8) then
!         write(6,*)'surfrdMod error: sum(pct) over maxpatch_pft is < 0.'
!         write(6,*)sumpct, nl
!         call endrun()
!       end if
!     end if

    end do   ! end of latitude loop

    ! Reset everithing if bare soil or special classes
    do nl = begg,endg
      if ( pctspecB(nl) > 99.0_r8 ) then
        pctpft_lunit(nl,:) = 0.0_r8
      end if
      if ( pctpft_lunit(nl,1) > 99.0_r8 ) then
        pctpft_lunit(nl,2:) = 0.0_r8
      end if
    end do

    ! Determine array [veg], which sets the PFT type for each of the [maxpatch]
    ! patches and array [wtxy], which sets the relative abundance of the PFT.
    ! Fill in PFTs for vegetated portion of grid cell. Fractional areas for
    ! these points [pctpft] pertain to "vegetated" area not to total grid area.
    ! So need to adjust them for fraction of grid that is vegetated.
    ! Next, fill in urban, lake, wetland, and glacier patches.

    do nl = begg,endg
      if (domain%pftm(nl) >= 0) then

        ! Naturally vegetated landunit

        do m = 1, maxpatch_pft
          vegxy(nl,m)  = pft(nl,m)
          wtxy(nl,m) = pctpft_lunit(nl,m) / 100.0_r8
#if (defined CN)
          ! the following test prevents the assignment of temperate deciduous
          ! vegetation types in the tropics
          ! 1. broadleaf deciduous temperate tree -> broadleaf deciduous tropical tree

          if (vegxy(nl,m) == 7 .and. abs(domain%latc(nl)) < 23.5_r8) vegxy(nl,m) = 6

          ! 2. broadleaf deciduous temperate shrub -> broadleaf deciduous tropical tree
          ! this reassignment from shrub to tree is necessary because there is currently no
          ! tropical deciduous broadleaf shrub type defined.

          if (vegxy(nl,m) == 10 .and. abs(domain%latc(nl)) < 23.5_r8) vegxy(nl,m) = 6
#endif
        end do
        ! Crop landunit

        if (create_crop_landunit) then
          do m = 1,maxpatch_cft
            vegxy(nl,npatch_glacier+m)  = cft(nl,m)
            wtxy(nl,npatch_glacier+m) = pctcft_lunit(nl,m) / 100._r8
          end do
        end if
      end if
    end do

    do nl = begg,endg
      if ( wtxy(nl,npatch_wet) > 0.99_r8 ) then
        wtxy(nl,:) = 0.0_r8
        wtxy(nl,npatch_wet) = 1.0_r8
      end if
      if ( wtxy(nl,npatch_glacier) > 0.99_r8 ) then
        wtxy(nl,:) = 0.0_r8
        wtxy(nl,npatch_glacier) = 1.0_r8
      end if
      if ( wtxy(nl,npatch_urban) > 0.99_r8 ) then
        wtxy(nl,:) = 0.0_r8
        wtxy(nl,npatch_urban) = 1.0_r8
      end if
      if ( wtxy(nl,npatch_lake) > 0.99_r8 ) then
        wtxy(nl,:) = 0.0_r8
        wtxy(nl,npatch_lake) = 1.0_r8
      end if
    end do

    found = .false.
    nindx(:) = -1
    sumvec(:) = 0._r8
    do nl = begg,endg
      if (domain%pftm(nl) >= 0) then
        sumvec(nl) = abs(sum(wtxy(nl,:)))
      end if
    end do
    do nl = begg,endg
      if (sumvec(nl) > 1.00005_r8 .and. domain%pftm(nl)>=0) then
        found = .true.
        nindx(nl) = nl
        exit
      endif
    end do
    if ( found ) then
      write (6,*)'surfrd error: WTXY > 1 occured'
      do nl = begg,endg
        if ( nindx(nl) > -1 ) then
          write (6,*) nl , sumvec(nl)
          do m = 1 , maxpatch
            write (6,*) m , wtxy(nl,m)
          end do
        end if
      end do
      call endrun()
    end if

    deallocate(sumvec,cft,pft)
    deallocate(pctcft_lunit,pctpft_lunit,pctpft)
    deallocate(wsti)
!abt below
    deallocate(adj_wst,adj_pctpft)
!abt above
  end subroutine surfrd_wtxy_veg_rank

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_veg_all
!
! !INTERFACE:
!  subroutine surfrd_wtxy_veg_all(ncid, pctspec, vegxy, wtxy, domain)
  subroutine surfrd_wtxy_veg_all(ncid, domain)
!
! !DESCRIPTION:
! Determine wtxy and veg arrays for non-dynamic landuse mode
!
! !USES:
    use domainMod   , only : domain_type
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer , intent(in)    :: ncid       ! netcdf file id 
!    real(r8), intent(in)    :: pctspec(:) ! percent wrt gcell of spec lunits
!    integer , intent(inout) :: vegxy(:,:)   ! PFT
!    real(r8), intent(inout) :: wtxy(:,:)  ! subgrid weights
    type(domain_type),intent(in) :: domain ! domain associated with wtxy
!
! !CALLED FROM:
! subroutine surfrd in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein, Sam Levis and Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m,mp7,mp8,mp11,n,nl,ns         ! indices
    integer  :: begg,endg                      ! beg/end gcell index
    integer  :: dimid,varid                    ! netCDF id's
    integer  :: start(4),count(4)              ! netcdf start/count arrays
    integer  :: ier                            ! error status	
    real(r8) :: sumpct                         ! sum of %pft over maxpatch_pft
    real(r8),allocatable :: pctpft(:,:)        ! percent of vegetated gridcell area for PFTs
    real(r8),pointer :: arrayl(:)              ! local array
    integer  :: ret, time_index
    real(r8) :: numpftp1data(0:numpft)         
    integer  :: strt3(3)                        ! Start index to read in
    integer  :: cnt3(3)                         ! Number of points to read in
    real(r8) :: closelat                        ! Single-column latitude value
    real(r8) :: closelon                        ! Single-column longitude value
    integer  :: closelatidx                     ! Single-column latitude index to retrieve
    integer  :: closelonidx                     ! Single-column longitude index to retrieve
    character(len=32) :: subname = 'surfrd_wtxy_veg_all'  ! subroutine name
!-----------------------------------------------------------------------

    ns = domain%ns
    call get_proc_bounds(begg,endg)
    allocate(pctpft(begg:endg,0:numpft))

    if (masterproc) then
!       call check_dim(ncid, 'lsmpft', numpft+1) abt
       call check_dim(ncid, 'level', numpft+1)
    endif

    allocate(arrayl(begg:endg))
    do n = 0,numpft
       start(1) = 1
       count(1) = domain%ni
       start(2) = 1
       count(2) = domain%nj
       start(3) = n+1         ! dataset is 1:numpft+1, not 0:numpft
       count(3) = 1

       start(4) = 1
       count(4) = 1     ! abt added to include time dimesion

       if (single_column) then
          call scam_setlatlonidx(ncid,scmlat,scmlon,closelat,closelon,closelatidx,closelonidx)
          start(1)=closelonidx
          start(2)=closelatidx
       endif
       call ncd_iolocal(ncid, 'PCT_PFT', 'read', arrayl, begg, endg, gsMap_lnd_gdc2glo, &
                        perm_lnd_gdc2glo, start, count)
       pctpft(begg:endg,n) = arrayl(begg:endg)
    enddo
    deallocate(arrayl)

    do nl = begg,endg
       if (domain%pftm(nl) >= 0) then

          ! Error check: make sure PFTs sum to 100% cover for vegetated landunit 
          ! (convert pctpft from percent with respect to gridcel to percent with 
          ! respect to vegetated landunit)

          if ( pctspecB(nl) < 100._r8 .and. pctspecB(nl) > 0.0 ) then
             sumpct = 0._r8
             do m = 0,numpft
                sumpct = sumpct + pctpft(nl,m)
             end do
             if (abs(sumpct - 100._r8 + pctspecB(nl) ) > 0.1e-8_r8) then
                write(6,*)'surfrdMod error: sum(pct) over numpft+1 is not = 100.'
                write(6,*) nl , sumpct, pctspecB(nl)
                call endrun()
             end if
             if (sumpct < -0.00000001_r8) then
                write(6,*)'surfrdMod error: sum(pct) over numpft+1 is < 0.'
                write(6,*) sumpct, nl
                call endrun()
             end if
          end if

          ! Set weight of each pft wrt gridcell (note that maxpatch_pft = numpft+1 here)

          do m = 1,numpft+1
             vegxy(nl,m)  = m-1
             wtxy(nl,m) = pctpft(nl,m) /100.0_r8
          end do

#if (defined CN)
          ! the following test prevents the assignment of temperate deciduous
          ! vegetation types in the tropics
          ! 1. broadleaf deciduous temperate tree (type7) -> broadleaf deciduous tropical tree (type6)
          ! N.B. the veg and wtxy arrays start at 1, so index 1 corresponds to
          ! veg type 0.  So in this case I want to trap on veg types 7 and 10, 
          ! which are indices 8 and 11. Moving to vegtype6, or index 7
          mp7 = 7
          mp8 = 8
          mp11 = 11
          if (abs(domain%latc(nl)) < 23.5_r8 .and. wtxy(nl,mp8) > 0._r8) then
             if (masterproc) then
                write(6,*)'surfrdMod warning: reassigning temperate tree -> tropical tree'
                write(6,*)'nl,lat,veg7wt,veg6wt,type'
                write(6,*) nl,domain%latc(nl),wtxy(nl,mp8),wtxy(nl,mp7),vegxy(nl,mp8)
             end if
             wtxy(nl,mp7) = wtxy(nl,mp7) + wtxy(nl,mp8)
             wtxy(nl,mp8) = 0._r8
             if (masterproc) then
                write(6,*) nl,domain%latc(nl),wtxy(nl,mp8),wtxy(nl,mp7)
             end if
          end if

          ! 2. broadleaf deciduous temperate shrub (type10) -> broadleaf deciduous tropical tree (type6)
          ! this reassignment from shrub to tree is necessary because there is currently no
          ! tropical deciduous broadleaf shrub type defined.

          if (abs(domain%latc(nl)) < 23.5_r8 .and. wtxy(nl,mp11) > 0._r8) then
             if (masterproc) then
                write(6,*)'surfrdMod warning: reassigning temperate shrub -> tropical tree'
                write(6,*)'nl,lat,veg10wt,veg6wt,type'
                write(6,*) nl,domain%latc(nl),wtxy(nl,mp11),wtxy(nl,mp7),vegxy(nl,mp11)
             end if
             wtxy(nl,mp7) = wtxy(nl,mp7) + wtxy(nl,mp11)
             wtxy(nl,mp11) = 0._r8
             if (masterproc) then
                write(6,*) nl,domain%latc(nl),wtxy(nl,mp11),wtxy(nl,mp7)
             end if
          end if
#endif
       end if
    end do

    deallocate(pctpft)

  end subroutine surfrd_wtxy_veg_all

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: surfrd_wtxy_veg_dgvm
!
! !INTERFACE:
!  subroutine surfrd_wtxy_veg_dgvm(pctspec, vegxy, wtxy, domain)
  subroutine surfrd_wtxy_veg_dgvm(domain)
!
! !DESCRIPTION:
! Determine wtxy and vegxy for DGVM mode.
!
! !USES:
    use pftvarcon   , only : crop, noveg
    use domainMod   , only : domain_type
!
! !ARGUMENTS:
    implicit none
!    real(r8), intent(in)    :: pctspec(:) ! percent gridcell of special landunits
!    integer , intent(inout) :: vegxy(:,:)   ! PFT
!    real(r8), intent(inout) :: wtxy(:,:)  ! subgrid weights
    type(domain_type),intent(in) :: domain ! domain associated with wtxy
!
! !CALLED FROM:
! subroutine surfrd in this module
!
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/04
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: m,nl         ! indices
    integer  :: begg,endg   ! beg/end gcell index
!-----------------------------------------------------------------------

    call get_proc_bounds(begg,endg)
    do nl = begg,endg
       do m = 1, maxpatch_pft
          vegxy(nl,m)  = noveg 
          wtxy(nl,m) = 1.0_r8/maxpatch_pft * (100._r8-pctspec(nl))/100._r8
       end do
    end do

  end subroutine surfrd_wtxy_veg_dgvm
   
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: surfrd_mkrank
!
! !INTERFACE:
  subroutine surfrd_mkrank (n, a, miss, iv, num)
!
! !DESCRIPTION:
! Return indices of largest [num] values in array [a]
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use abortutils, only : endrun
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: n        ! array length
    real(r8), intent(in) :: a(0:n)   ! array to be ranked
    integer , intent(in) :: miss     ! missing data value
    integer , intent(in) :: num      ! number of largest values requested
    integer , intent(out):: iv(num)  ! index to [num] largest values in array [a]
!
! !CALLED FROM:
! ! subroutine surfrd_wtxy_veg_rank in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
!
!EOP
!
! !LOCAL VARIABLES:
    real(r8) :: a_max       ! maximum value in array
    integer  :: i           ! array index
    real(r8) :: delmax      ! tolerance for finding if larger value
    integer  :: m           ! do loop index
    integer  :: k           ! do loop index
    logical  :: exclude     ! true if data value has already been chosen
!-----------------------------------------------------------------------

    delmax = 1.e-06_r8

    ! Find index of largest non-zero number
    
    iv(1) = miss
    a_max = -9999._r8

    do i = 0, n
       if (a(i)>0._r8 .and. (a(i)-a_max)>delmax) then
          a_max = a(i)
          iv(1)  = i
       end if
    end do

    ! iv(1) = miss indicates no values > 0. this is an error

    if (iv(1) == miss) then
       write (6,*) 'surfrd_mkrank error: iv(1) = missing'
       call endrun
    end if

    ! Find indices of the next [num]-1 largest non-zero number.
    ! iv(m) = miss if there are no more values > 0

    do m = 2, num
       iv(m) = miss
       a_max = -9999._r8
       do i = 0, n

          ! exclude if data value has already been chosen

          exclude = .false.
          do k = 1, m-1
             if (i == iv(k)) exclude = .true.
          end do

          ! if not already chosen, see if it is the largest of
          ! the remaining values

          if (.not. exclude) then
             if (a(i)>0._r8 .and. (a(i)-a_max)>delmax) then
                a_max = a(i)
                iv(m)  = i
             end if
          end if
       end do
    end do

  end subroutine surfrd_mkrank


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pft_adjustment
!
! !INTERFACE:
!  subroutine pft_adjustment(begg,endg,numpft,pftpct,domain)
  subroutine pft_adjustment(begg,endg,ipft,nns,pft_wst,pft_wst_sum,pft_pctpft,domain)
!
! !DESCRIPTION:
! Checks to see if there are gridcells with positive pftmask but no
! percent pft cover.  If there is no vegetation cover in a gridcell
! with a pftmask = 1 then use the BATS landuse type and assign the
! equivalent CLM pft type (using the technique in Bonan et al 2002)
!
! !USES:
    use domainMod   , only : domain_type 
    use decompMod   , only : adecomp
    use clm_varsur  , only : wtxy,satbrt_clm
    use mod_clm
    use mod_dynparam
!
    implicit none
    include 'netcdf.inc'
!
! !ARGUMENTS:
    integer  , intent(in)     :: nns                   !
    integer  , intent(in)     :: endg                  !
    integer  , intent(in)     :: begg                  !
    integer  , intent(out)    :: ipft                  !
    real(r8) , intent(out)    :: pft_wst_sum           !
    real(r8) , intent(out)    :: pft_wst(0:numpft)            ! total number of pfts 
    real(r8) , intent(out)    :: pft_pctpft(begg:endg,0:numpft)       ! percent pft cover 
    type(domain_type),intent(in) :: domain               ! domain type structure
!
! !CALLED FROM:
! subroutine surfrd_wtxy_veg_rank in this module
!
! !REVISION HISTORY:
! Created by Ahmed Tawfik
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: nni,ni,nj                       ! size of grid on file
    integer :: ncid,dimid,varid                ! netCDF id's
    integer :: nnj                             ! error status
    integer :: n                               ! counters
    real(r8) :: summ                           ! Number of points to read in
    real(r8),allocatable :: pft2d(:,:)         ! temp array for pftmask calc
    real(r8),allocatable :: sumpft(:)          ! landfraction read in from RCMnavy abt
    real(r8) :: bats_lu                        ! temporary array from BATS landuse
    character(len=32) :: subname = 'pft_adjustment'     ! subroutine name
!-----------------------------------------------------------------------

     pft_pctpft(nns,:) = 0._r8
     pft_wst(:)        = 0._r8          
     ipft = 0
          
     nni = adecomp%gdc2i(nns)
     nnj = adecomp%gdc2j(nns)
     if(satbrt_clm(nni,nnj).lt.0.9 .or. satbrt_clm(nni,nnj).gt.22.1) then
       write(*,*)"BATS landuse correction is not set properly :::: at nl = ",nns
       write(*,*)"SATBRT   = ",satbrt_clm(nni,nnj)
       write(*,*)"lat/lat  = ",nni,nnj
       call endrun()
     else
       bats_lu = satbrt_clm(nni,nnj)
     endif
       
     pft_wst_sum = 100.0_r8 - pctspecB(nns)

     ! Desert/Tundra/semidesert from BATS equals no veg for CLM
     if(bats_lu.eq.8 .or. bats_lu.eq.9 .or. bats_lu.eq.11) then
       pft_pctpft(nns,0) = pft_wst_sum
       pft_wst(0)        = pft_wst_sum
     endif
     
     ! Short/Tall grass from BATS equals grass for CLM
     if(bats_lu.eq.2 .or. bats_lu.eq.7) then
       pft_pctpft(nns,14) = pft_wst_sum
       pft_wst(14)        = pft_wst_sum
     endif

     ! Evergreen Needleleaf from BATS eq NeEvTemp and Boreal
     if(bats_lu.eq.3) then
       if(domain%latc(nns).gt.47) then  !boreal
         pft_pctpft(nns,2) = pft_wst_sum
         pft_wst(2)        = pft_wst_sum
       else                             !temperate
         pft_pctpft(nns,1) = pft_wst_sum
         pft_wst(1)        = pft_wst_sum
       endif
     endif
    
     ! Decidious Needleleaf from BATS eq NeDeBo
     if(bats_lu.eq.4) then
       pft_pctpft(nns,3) = pft_wst_sum
       pft_wst(3)        = pft_wst_sum
     endif 

     ! Deciduous Broadleaf from BATS eq Broadleaf Deciduous tropical/temperate/boreal
     if(bats_lu.eq.5) then
       if(domain%latc(nns).gt.47) then  !boreal
         pft_pctpft(nns,8) = pft_wst_sum
         pft_wst(8)        = pft_wst_sum
       elseif(domain%latc(nns).lt.15 .and. &
              domain%latc(nns).gt.-15) then     !tropical
         pft_pctpft(nns,6) = pft_wst_sum
         pft_wst(6)        = pft_wst_sum
       else                              !temperate
         pft_pctpft(nns,7)  = pft_wst_sum
         pft_wst(7)         = pft_wst_sum
       endif
     endif

     ! Evergreen Broadleaf from BATS eq Broadleaf evergreen tropical/temperate
     if(bats_lu.eq.6) then
       if(domain%latc(nns).lt.15 .and. &
          domain%latc(nns).gt.-15) then     !tropical
         pft_pctpft(nns,4) = pft_wst_sum
         pft_wst(4)        = pft_wst_sum
       else                              !temperate
         pft_pctpft(nns,5)  = pft_wst_sum
         pft_wst(5)         = pft_wst_sum
       endif
     endif
       
     ! Evergreen Shrub from BATS eq Broadleaf Evergreen Shrub
     if(bats_lu.eq.16) then
       pft_pctpft(nns,9) = pft_wst_sum
       pft_wst(9)        = pft_wst_sum
     endif

     ! Deciduous Shrub from BATS eq Broadleaf Deciduous temperate/boreal shurbs
     if(bats_lu.eq.17) then
       if(domain%latc(nns).gt.47) then  !boreal
         pft_pctpft(nns,11)  = pft_wst_sum
         pft_wst(11)         = pft_wst_sum
       else                              !temperate
         pft_pctpft(nns,10)  = pft_wst_sum
         pft_wst(10)         = pft_wst_sum
       endif
     endif

     ! Mixed forest from BATS eq some percent of needleleaf and broadleaf
     if(bats_lu.eq.18) then
       if(domain%latc(nns).gt.47) then  !boreal
         pft_pctpft(nns,8)  = pft_wst_sum/2.0_r8
         pft_pctpft(nns,2)  = pft_pctpft(nns,8)/2.0_r8
         pft_pctpft(nns,3)  = pft_wst_sum - &
           (pft_pctpft(nns,8) + pft_pctpft(nns,2))
         pft_wst(8)         = pft_pctpft(nns,8)
         pft_wst(2)         = pft_pctpft(nns,2)
         pft_wst(3)         = pft_pctpft(nns,3)
       elseif(domain%latc(nns).lt.15 .and. &
              domain%latc(nns).gt.-15) then     !tropical
         pft_pctpft(nns,4)  = pft_wst_sum/2.0_r8
         pft_pctpft(nns,6)  = pft_wst_sum - pft_pctpft(nns,4)
         pft_wst(4)         = pft_pctpft(nns,4)
         pft_wst(6)         = pft_pctpft(nns,6)
       else                              !temperate
         pft_pctpft(nns,1)  = pft_wst_sum/3.0_r8
         pft_pctpft(nns,5)  = pft_pctpft(nns,1)
         pft_pctpft(nns,7)  = pft_wst_sum - &
                       (pft_pctpft(nns,5)+pft_pctpft(nns,1))
         pft_wst(7)         = pft_pctpft(nns,7)
         pft_wst(5)         = pft_pctpft(nns,5)
         pft_wst(1)         = pft_pctpft(nns,1)
       endif
     endif

     ! Forest mosaic from BATS eq some forests and grass
     if(bats_lu.eq.19) then
       if(domain%latc(nns).gt.47) then  !boreal
         pft_pctpft(nns,8)  = pft_wst_sum/6.0_r8
         pft_pctpft(nns,2)  = pft_pctpft(nns,8)
         pft_pctpft(nns,3)  = pft_pctpft(nns,2)
         pft_pctpft(nns,14) = pft_wst_sum - &
           (pft_pctpft(nns,8)+pft_pctpft(nns,2)+pft_pctpft(nns,3))
         pft_wst(8)         = pft_pctpft(nns,8)
         pft_wst(2)         = pft_pctpft(nns,2)
         pft_wst(3)         = pft_pctpft(nns,3)
         pft_wst(14)        = pft_pctpft(nns,14)
       elseif(domain%latc(nns).lt.15 .and. domain%latc(nns).gt.-15) then     !tropical
         pft_pctpft(nns,4)  = pft_wst_sum/4.0_r8
         pft_pctpft(nns,6)  = pft_pctpft(nns,4)
         pft_pctpft(nns,14) = pft_wst_sum - &
            (pft_pctpft(nns,4) + pft_pctpft(nns,6))
         pft_wst(4)         = pft_pctpft(nns,4)
         pft_wst(6)         = pft_pctpft(nns,6)
         pft_wst(14)        = pft_pctpft(nns,14)
       else                              !temperate
         pft_pctpft(nns,1)  = pft_wst_sum/6.0_r8
         pft_pctpft(nns,5)  = pft_pctpft(nns,1)
         pft_pctpft(nns,7)  = pft_pctpft(nns,5)
         pft_pctpft(nns,14) = pft_wst_sum - &
           (pft_pctpft(nns,1) + pft_pctpft(nns,5) + pft_pctpft(nns,7))
         pft_wst(1)         = pft_pctpft(nns,1)
         pft_wst(5)         = pft_pctpft(nns,5)
         pft_wst(7)         = pft_pctpft(nns,7)
         pft_wst(14)        = pft_pctpft(nns,14)
      endif
    endif

    ! Crop/irrigated from BATS eq half corn and wheat mix
    if(bats_lu.eq.1 .or. bats_lu.eq.10) then
      pft_pctpft(nns,15)  = pft_wst_sum/2.0_r8
      pft_pctpft(nns,16)  = pft_wst_sum - pft_pctpft(nns,15)
      pft_wst(15)         = pft_pctpft(nns,15)
      pft_wst(16)         = pft_pctpft(nns,16)
    endif           

    ! Lake from BATS equals lake for CLM
    if(bats_lu.eq.14) then
      pctspec(nns)       = 100._r8
      pctspecB(nns)      = 100._r8
      wtxy(nns,npatch_lake) = 1.0_r8  
      wtxy(nns,npatch_glacier) = 0.0_r8  
      wtxy(nns,npatch_wet) = 0.0_r8  
      wtxy(nns,npatch_urban) = 0.0_r8  
    endif

    ! Glacier from BATS equals glacier for CLM
    if(bats_lu.eq.12) then
      pft_wst_sum        = 0._r8
      pctspec(nns)       = 100._r8
      pctspecB(nns)      = 100._r8
      wtxy(nns,npatch_glacier) = 1.0_r8  
      wtxy(nns,npatch_lake) = 0.0_r8  
      wtxy(nns,npatch_wet) = 0.0_r8  
      wtxy(nns,npatch_urban) = 0.0_r8  
    endif

    ! Wetland from BATS equals wetland for CLM
    if(bats_lu.eq.13) then
      pctspec(nns)       = 100._r8
      pctspecB(nns)      = 100._r8
      wtxy(nns,npatch_wet) = 1.0_r8  
      wtxy(nns,npatch_glacier) = 0.0_r8  
      wtxy(nns,npatch_lake) = 0.0_r8  
      wtxy(nns,npatch_urban) = 0.0_r8  
    endif

    ! Urban/Semi Urban from BATS equals urban for CLM
    if(bats_lu.gt.20) then
      pctspec(nns)       = 100._r8
      pctspecB(nns)      = 100._r8
      wtxy(nns,npatch_urban) = 1.0_r8  
      wtxy(nns,npatch_wet) = 0.0_r8  
      wtxy(nns,npatch_glacier) = 0.0_r8  
      wtxy(nns,npatch_lake) = 0.0_r8  
    endif

    ! Bare ground for any other landtype 
    if(bats_lu.eq.15 .or. bats_lu.eq.20) then
      pft_wst_sum        = 100._r8
      pft_pctpft(nns,0)  = 100._r8
      pft_wst(0)         = 100._r8
      wtxy(nns,npatch_wet)     = 0.0_r8
      wtxy(nns,npatch_urban)   = 0.0_r8
      wtxy(nns,npatch_lake)    = 0.0_r8
      wtxy(nns,npatch_glacier) = 0.0_r8
    endif

  end subroutine pft_adjustment
!
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm2bats_conversion
!
! !INTERFACE:
  subroutine clm2bats_conversion(ncidpft,ncidlak,ncidglac,ncidurb)
!
! !DESCRIPTION:
! 1) Uses the maximum land cover type in CLM to describe the land cover 
! type in terms of BATS land use categories.  In some cases use some
! percent threshold limits to make the conversion from CLM land use to
! to BATS land use.  For example, BATS cropland land use type would be
! equivalent to having the sum of wheat and corn being greater than
! 50%.  
!
! Reference for vegetation cover conversion:
! Bonan GB, Levis S,Kergoat L et al.(2002a) Landscapes as
! patches of plant functional types: an integrating concept for 
! climate and ecosystem models. Global Biogeochemical Cycles, 
! 16, 1021
!
! !USES:
    use clm_varsur  , only : clm2bats_veg,clm_fracveg,satbrt_clm
    use clm_varsur  , only : landmask
    use mod_clm
    use mod_dynparam
!
    implicit none
    include 'netcdf.inc'
!
! !ARGUMENTS:
    integer , intent(in)    :: ncidlak      ! netcdf file id for lake/wetland
    integer , intent(in)    :: ncidglac     ! netcdf file id for glacier
    integer , intent(in)    :: ncidurb      ! netcdf file id for urban
    integer , intent(in)    :: ncidpft      ! netcdf file id for soil
!
! !CALLED FROM:
! subroutine rcmsurf in this module
!
! !REVISION HISTORY:
! Created by Ahmed Tawfik
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: nni,ni                         ! longitude index
    integer  :: nnj,nj                         ! latitude index
    integer  :: n,nns,ier                      ! counters
    integer  :: varid,ncid,strt3(4),cnt3(4)    !
    real(r8) :: summ                           ! temporary
    real(r8) :: sum2                           ! temporary 
    real(r8) :: pctveg                         ! percent vegetation in gridcell
    real(r8) :: purb                           ! percent urban in gridcell
    real(r8) :: plake                          ! percent lake in gridcell
    real(r8) :: pwet                           ! percent wetland in gridcell
    real(r8) :: pgla                           ! percent glacier in gridcell
    real(r8), allocatable :: perpft(:,:)
    real(r8), allocatable :: pcturb(:)
    real(r8), allocatable :: pctlak(:)
    real(r8), allocatable :: pctwet(:)
    real(r8), allocatable :: pctgla(:)
    real(r8), allocatable :: pctpft(:)
    real(r8), allocatable :: array1(:)
    character(len=32) :: subname = 'clm2bats_conversion'     ! subroutine name
!-----------------------------------------------------------------------

    if(masterproc) then

!******* read in landunit variables
      allocate(pcturb(lsmlat*lsmlon),pctwet(lsmlat*lsmlon))  
      allocate(pctgla(lsmlat*lsmlon),pctlak(lsmlat*lsmlon))  
      allocate(perpft(lsmlat*lsmlon,0:numpft))
      allocate(pctpft(0:numpft))

      strt3(1) = 1
      cnt3(1)  = lsmlon
      strt3(2) = 1
      cnt3(2)  = lsmlat
      strt3(3) = 1
      cnt3(3)  = 1
      cnt3(4)  = 1
      strt3(4) = 1   !abt added this to include time dimension

      ncid = ncidurb
      call check_ret(nf_inq_varid(ncid, 'PCT_URBAN', varid), subname)
      call check_ret(nf_get_vara_double(ncid, varid, strt3(1:3), cnt3(1:3), pcturb), subname)

      ncid = ncidlak
      call check_ret(nf_inq_varid(ncid, 'PCT_WETLAND', varid), subname)
      call check_ret(nf_get_vara_double(ncid, varid, strt3(1:3), cnt3(1:3), pctwet), subname)

      call check_ret(nf_inq_varid(ncid, 'PCT_LAKE', varid), subname)
      call check_ret(nf_get_vara_double(ncid, varid, strt3(1:3), cnt3(1:3), pctlak), subname)

      ncid = ncidglac 
      call check_ret(nf_inq_varid(ncid, 'PCT_GLACIER', varid), subname)
      call check_ret(nf_get_vara_double(ncid, varid, strt3(1:3), cnt3(1:3), pctgla), subname)

      ncid = ncidpft
      allocate(array1(lsmlat*lsmlon))
      do n = 0,numpft
        strt3(3) = n+1 ! dataset is 1:numpft+1, not 0:numpft
        cnt3(3)  = 1
        cnt3(4)  = 1
        strt3(4) = 1   !abt added this to include time dimension

        call check_ret(nf_inq_varid(ncid, 'PCT_PFT', varid), subname)
        call check_ret(nf_get_vara_double(ncid, varid, strt3, cnt3, array1), subname)

        perpft(:,n) = array1(:)
      enddo
      deallocate(array1)


!******* Figure out what land use type covers the most area per gridcell
   
      do nnj = 1,lsmlat
      do nni = 1,lsmlon
         nns = (nnj-1)*lsmlon + nni


       ! Make sure only land points are counted
        if(landmask(nni,nnj).lt.0.0001) then
           clm2bats_veg(nni,nnj) = 0
           clm_fracveg(nni,nnj)  = 0

        else
           purb  = pcturb(nns)
           plake = pctlak(nns)
           pgla  = pctgla(nns)
           pwet  = pctwet(nns)

           pctveg = 100._r8 - (purb + pwet + plake + pgla)
           pctpft(:) = perpft(nns,:) 

           clm_fracveg(nni,nnj) = (((100._r8 - pctpft(0))/100._r8)*pctveg)/100._r8
           clm2bats_veg(nni,nnj) = satbrt_clm(nni,nnj)   

 
         ! clm2bats_veg is initialized to all ocean and the respective 
         ! land description is filled in below. clm2bats_veg is
         ! initialized init_clm(ser/para).F

         ! Urban land in CLM equates to nothing and is currently not
         ! parameterized so set to BATS irrigated crop for now
           if( purb.gt.pctveg .and. purb.gt.pwet .and. purb.gt.plake &
               .and. purb.gt.pgla) then
               
               clm2bats_veg(nni,nnj) = 10
           end if            



         ! Wetland in CLM equates to BATS bogs/marsh type (number 13)
           if( pwet.gt.pctveg .and. pwet.gt.purb .and. pwet.gt.plake &
               .and. pwet.gt.pgla) then
               
               clm2bats_veg(nni,nnj) = 13
           end if            



         ! Glacier in CLM equates to BATS glacier/ice cap type (number 12)
           if( pgla.gt.pctveg .and. pgla.gt.purb .and. pgla.gt.plake &
               .and. pgla.gt.pwet) then
               
               clm2bats_veg(nni,nnj) = 12
           end if            



         ! Lake in CLM equates to BATS inland water type (number 14)
           if( plake.gt.pctveg .and. plake.gt.purb .and. plake.gt.pwet &
               .and. plake.gt.pgla) then
               
               clm2bats_veg(nni,nnj) = 14
           end if            


         ! Vegetation is the dominant land cover type
         ! If it is the most dominant then we must discriminate
         ! to only use the correct BATS type

           if( pctveg.gt.pwet .and. pctveg.gt.purb .and. pctveg.gt.plake &
               .and. pctveg.gt.pgla) then

              ! Sum of wheat and corn from CLM equates to BATS cropland  
               summ = pctpft(15) + pctpft(16)
               sum2 = sum(pctpft(:)) - summ
               if(summ .gt. sum2) then
                 clm2bats_veg(nni,nnj)  = 1
               end if


              ! Non-artic grass from CLM equates to BATS short grass  
               summ = pctpft(13)
               sum2 = sum(pctpft(:)) - summ
               if(summ .gt. sum2) then
                 clm2bats_veg(nni,nnj)  = 2
               end if


              ! Sum of all needleleaf evergreens from CLM equates to
              ! same in BATS   
               summ = pctpft(1) + pctpft(2)
               sum2 = sum(pctpft(:)) - summ
               if(summ .gt. sum2) then
                 clm2bats_veg(nni,nnj)  = 3
               end if


              ! Deciduous needleleaf from CLM equates to BATS
              ! deciduous needleleaf  
               summ = pctpft(3)
               sum2 = sum(pctpft(:)) - summ
               if(summ .gt. sum2) then
                 clm2bats_veg(nni,nnj)  = 4
               end if


              ! Deciduous broadleaf from CLM equates to BATS
              ! deciduous broadleaf  
               summ = pctpft(6) + pctpft(7) + pctpft(8)
               sum2 = sum(pctpft(:)) - summ
               if(summ .gt. sum2) then
                 clm2bats_veg(nni,nnj)  = 5
               end if

              ! Evergreen Broadleaf from CLM equates to BATS
              ! evergreen broadleaf  
               summ = pctpft(4) + pctpft(5)
               sum2 = sum(pctpft(:)) - summ
               if(summ .gt. sum2) then
                 clm2bats_veg(nni,nnj)  = 6
               end if


              ! C4 grasses from CLM equates to BATS
              ! Tall grass  
               summ = pctpft(15)
               sum2 = sum(pctpft(:)) - summ
               if(summ .gt. sum2) then
                 clm2bats_veg(nni,nnj)  = 7
               end if


              ! Bare Ground from CLM equates to BATS
              ! Desert  
               summ = pctpft(0)
               sum2 = pgla + plake + pwet + purb
               if(summ.ge.98 .and. sum2.lt.95) then
                 clm2bats_veg(nni,nnj)  = 8
               end if


              ! Artic Grass from CLM equates to BATS
              ! Tundra  
               summ = pctpft(12)
               sum2 = sum(pctpft(:)) - summ
               if(summ .gt. sum2) then
                 clm2bats_veg(nni,nnj)  = 9
               end if


              ! Bare ground percent from CLM equates to BATS
              ! semi-desert  
               summ = pctpft(0)
               sum2 = sum(pctpft(:)) - summ
               if(summ.gt.90 .and. summ.lt.98 ) then
                 clm2bats_veg(nni,nnj)  = 11
               end if


              ! Broadleaf Evergreen Shrub from CLM equates to BATS
              ! Evergreen Shrub  
               summ = pctpft(9)
               sum2 = sum(pctpft(:)) - summ
               if(summ .gt. sum2) then
                 clm2bats_veg(nni,nnj)  = 16
               end if


              ! Sum of Deciduous Shrub from CLM equates to BATS
              ! Deciduous Shrub  
               summ = pctpft(10) + pctpft(11)
               sum2 = sum(pctpft(:)) - summ
               if(summ .gt. sum2) then
                 clm2bats_veg(nni,nnj)  = 17
               end if
               
           end if      ! vegetation dominant gridcell if statement

        end if  ! landmask if statement
      
      end do
      end do

      deallocate(pcturb,pctlak,pctgla,pctwet)
      deallocate(perpft)

    end if !masterproc

    call mpi_bcast (clm2bats_veg, size(clm2bats_veg), MPI_REAL8, 0, mpicom, ier)
    call mpi_bcast (clm_fracveg, size(clm_fracveg), MPI_REAL8, 0, mpicom, ier)



  end subroutine clm2bats_conversion



!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: clm_getsoitex
!
! !INTERFACE:
  subroutine clm_getsoitex()
!
! !DESCRIPTION:
! Subrountine only called when RegCM Aerosol/Dust module is on
! Used only for the RegCM Dust scheme (Zakey Dust Scheme)
! Divide the soil texture fractions into soil texture classes
! defined in inidust.F.  
!
! !USES:
    use domainMod   , only : domain_type 
    use clm_varsur  , only : clm_soitex
    use clm_varctl  , only : mksrf_fsoitex
    use fileutils   , only : getfil
    use mod_clm
    use mod_dynparam
!
    implicit none
    include 'netcdf.inc'
!
! !ARGUMENTS:
!
! !CALLED FROM:
! subroutine rcmsurf in this module
!
! !REVISION HISTORY:
! Created by Ahmed Tawfik
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: nni,nnj,ns,n,ni,nj
    integer  :: nlat,nlon
    integer  :: varid,ncid
    integer  :: start(4),count(4)
    real(r8) :: slt                                    ! temporary variable for Silt
    logical  :: vpresent
    real(r8), allocatable :: sand1d(:)                 ! used for percent sand
    real(r8), allocatable :: clay1d(:)                 ! used for percent clay
    character(len=256) :: locfn                        ! local filename
    character(len=32) :: subname = 'clm_getsoitex'     ! subroutine name
!-----------------------------------------------------------------------

            nlat = iy
            nlon = jx

            start(1) = 1
            start(2) = 1
            start(3) = 1
            start(4) = 1
            count(1) = nlon
            count(2) = nlat
            count(3) = 1
            count(4) = 1
            
            call getfil ("../Input/RCMsoitex.10level.nc", locfn, 0)
            call check_ret(nf_open(locfn, 0, ncid), subname)

            if(.not.allocated(clm_soitex)) allocate(clm_soitex(nlat,nlon))
            if(.not.allocated(sand1d)) allocate(sand1d(nlat*nlon))
            if(.not.allocated(clay1d)) allocate(clay1d(nlat*nlon))

            do n = 1,1
              start(3) = n
              call check_ret(nf_inq_varid(ncid, 'PCT_SAND', varid), subname)
              call check_ret(nf_get_vara_double(ncid,varid, &
              start,count,sand1d), subname)

              call check_ret(nf_inq_varid(ncid, 'PCT_CLAY', varid), subname)
              call check_ret(nf_get_vara_double(ncid,varid, &
              start,count,clay1d), subname)
            enddo

            do nj = 1,nlat
            do ni = 1,nlon
               ns = (nj-1)*nlon + ni

              slt = 100._r8 - (sand1d(ns) + clay1d(ns)) 

             ! Class 1
              if(sand1d(ns).ge.90) then
                 clm_soitex(nj,ni) = 1                  

 
             ! Classes 2 3 11 or 10
              elseif(sand1d(ns).ge.60 .and. sand1d(ns).lt.89) then
                 if(slt.lt.2) then
                   if(clay1d(ns).lt.35) clm_soitex(nj,ni) = 10 
                   if(clay1d(ns).ge.35) clm_soitex(nj,ni) = 11  
                 else
                   if(clay1d(ns).lt.5)  clm_soitex(nj,ni) = 2 
                   if(clay1d(ns).ge.5 .and. clay1d(ns).le.10) clm_soitex(nj,ni) = 3 
                 end if


             ! Classes 4 and 12
              elseif(sand1d(ns).ge.50 .and. sand1d(ns).lt.60) then
                 if(slt.lt.2) then
                   clm_soitex(nj,ni) = 12
                 else
                   clm_soitex(nj,ni) = 4  
                 end if


             ! Classes 5
              elseif(sand1d(ns).ge.45 .and. sand1d(ns).lt.50) then
                clm_soitex(nj,ni) = 5    


             ! Classes 6
              elseif(sand1d(ns).ge.35 .and. sand1d(ns).lt.45) then
                clm_soitex(nj,ni) = 6    


             ! Classes 7 and 8
              elseif(sand1d(ns).ge.30 .and. sand1d(ns).lt.35) then
                 if(clay1d(ns).gt.15 .and. clay1d(ns).lt.23) then
                   clm_soitex(nj,ni) = 7   
                 elseif(clay1d(ns).ge.23 .and. clay1d(ns).lt.30) then
                   clm_soitex(nj,ni) = 8 
                 end if


             ! Classes 6
              elseif(sand1d(ns).ge.0.0000000000005 .and. sand1d(ns).lt.30) then
                clm_soitex(nj,ni) = 9    


              end if
              

            end do ! end of gridcell loop 
            end do ! end of gridcell loop 

            if(allocated(sand1d)) deallocate(sand1d)
            if(allocated(clay1d)) deallocate(clay1d)
            call check_ret(nf_close(ncid), subname)


     end subroutine clm_getsoitex


end module surfrdMod
