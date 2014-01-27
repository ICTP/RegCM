module mod_clm_domain
  !
  ! Module containing 2-d global surface boundary data information
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mppparam
  use mod_runparams
  use mod_mpmessage

  implicit none

  private

  public :: domain_type

  !--- this typically contains local domain info with arrays dim begg:endg ---
  type domain_type
    integer(ik4) :: ns         ! global size of domain
    integer(ik4) :: ni , nj    ! global axis if 2d (nj=1 if unstructured)
    logical :: isgrid2d   ! true => global grid is lat/lon
    integer(ik4) :: nbeg , nend  ! local beg/end indices
    character(len=8) :: clmlevel   ! grid type
    ! land mask: 1 = land, 0 = ocean
    integer(ik4) , pointer , dimension(:) :: mask
    ! fractional land
    real(rk8) , pointer , dimension(:) :: frac
    ! topography
    real(rk8) , pointer , dimension(:) :: topo
    ! latitude of grid cell (deg)
    real(rk8) , pointer , dimension(:) :: latc
    ! longitude of grid cell (deg)
    real(rk8) , pointer , dimension(:) :: lonc
    ! grid cell area (km**2)
    real(rk8) , pointer , dimension(:) :: area
    ! pft mask: 1=real, 0=fake, -1=notset
    integer(ik4) , pointer , dimension(:) :: pftm
    ! glc mask: 1=sfc mass balance required by GLC component
    ! 0=SMB not required (default)
    integer(ik4) , pointer , dimension(:) :: glcmask
    character(len=16) :: set ! flag to check if domain is set
    logical :: decomped   ! decomposed locally or global copy
  end type domain_type

  type(domain_type) , public :: ldomain
  ! 1d lat/lons for 2d grids
  real(rk8) , allocatable , public , dimension(:) :: lon1d
  real(rk8) , allocatable , public , dimension(:) :: lat1d

  public :: domain_init          ! allocates/nans domain types
  public :: domain_clean         ! deallocates domain types
  public :: domain_check         ! write out domain info

  character(len=16) , parameter :: set   = 'domain_set      '
  character(len=16) , parameter :: unset = 'NOdomain_unsetNO'

  contains
  !
  ! This subroutine allocates and nans the domain type
  !
  subroutine domain_init(domain,isgrid2d,ni,nj,nbeg,nend,clmlevel)
    implicit none
    type(domain_type) :: domain           ! domain datatype
    logical , intent(in) :: isgrid2d      ! true => global grid is lat/lon
    integer(ik4) , intent(in) :: ni , nj  ! grid size, 2d
    integer(ik4) , intent(in) , optional :: nbeg , nend  ! beg/end indices
    character(len=*) , intent(in) , optional  :: clmlevel   ! grid type
    integer(ik4) :: ier
    integer(ik4) :: nb , ne

    nb = 1
    ne = ni*nj
    if ( present(nbeg) ) then
      if ( present(nend) ) then
        nb = nbeg
        ne = nend
      end if
    end if

    if ( domain%set == set ) then
      call domain_clean(domain)
    end if
    allocate(domain%mask(nb:ne),domain%frac(nb:ne),domain%latc(nb:ne), &
             domain%pftm(nb:ne),domain%area(nb:ne),domain%lonc(nb:ne), &
             domain%topo(nb:ne),domain%glcmask(nb:ne),stat=ier)
    if (ier /= 0) then
      call fatal(__FILE__,__LINE__, &
         'domain_init ERROR: allocate mask, frac, lat, lon, area ')
    end if

    if ( present(clmlevel) ) then
      domain%clmlevel = clmlevel
    end if

    domain%isgrid2d = isgrid2d
    domain%ns       = ni*nj
    domain%ni       = ni
    domain%nj       = nj
    domain%nbeg     = nb
    domain%nend     = ne
    domain%mask     = -9999
    domain%frac     = -1.0e36
    domain%topo     = 0.0D0
    domain%latc     = nan
    domain%lonc     = nan
    domain%area     = nan

    domain%set      = set
    if ( domain%nbeg == 1 .and. domain%nend == domain%ns ) then
      domain%decomped = .false.
    else
      domain%decomped = .true.
    end if

    domain%pftm     = -9999
    domain%glcmask  = 0  
  end subroutine domain_init
  !
  ! This subroutine deallocates the domain type
  !
  subroutine domain_clean(domain)
    implicit none
    type(domain_type) :: domain        ! domain datatype
    integer(ik4) :: ier
    if ( domain%set == set ) then
      if (myid == italk) then
        write(stdout,*) 'domain_clean: cleaning ',domain%ni,domain%nj
      end if
      deallocate(domain%mask,domain%frac,domain%latc, &
                 domain%lonc,domain%area,domain%pftm, &
                 domain%topo,domain%glcmask,stat=ier)
      if (ier /= 0) then
        call fatal(__FILE__,__LINE__, &
            'domain_clean ERROR: deallocate mask, frac, lat, lon, area ')
      end if
    else
      if (myid == italk) then
        write(stdout,*) 'domain_clean WARN: clean domain unecessary '
      end if
    end if
    domain%clmlevel   = unset
    domain%ns         = bigint
    domain%ni         = bigint
    domain%nj         = bigint
    domain%nbeg       = bigint
    domain%nend       = bigint
    domain%set        = unset
    domain%decomped   = .true.
  end subroutine domain_clean
  !
  ! This subroutine write domain info
  !
  subroutine domain_check(domain)
    implicit none
    type(domain_type) , intent(in) :: domain        ! domain datatype
    if ( myid == italk ) then
      write(stdout,*) '  domain_check set       = ',trim(domain%set)
      write(stdout,*) '  domain_check decomped  = ',domain%decomped
      write(stdout,*) '  domain_check ns        = ',domain%ns
      write(stdout,*) '  domain_check ni,nj     = ',domain%ni,domain%nj
      write(stdout,*) '  domain_check clmlevel  = ',trim(domain%clmlevel)
      write(stdout,*) '  domain_check nbeg,nend = ',domain%nbeg,domain%nend
      write(stdout,*) '  domain_check lonc      = ', &
              minval(domain%lonc),maxval(domain%lonc)
      write(stdout,*) '  domain_check latc      = ', &
              minval(domain%latc),maxval(domain%latc)
      write(stdout,*) '  domain_check mask      = ', &
              minval(domain%mask),maxval(domain%mask)
      write(stdout,*) '  domain_check frac      = ', &
              minval(domain%frac),maxval(domain%frac)
      write(stdout,*) '  domain_check topo      = ', &
              minval(domain%topo),maxval(domain%topo)
      write(stdout,*) '  domain_check area      = ', &
              minval(domain%area),maxval(domain%area)
      write(stdout,*) '  domain_check pftm      = ', &
              minval(domain%pftm),maxval(domain%pftm)
      write(stdout,*) '  domain_check glcmask   = ', &
              minval(domain%glcmask),maxval(domain%glcmask)
      write(stdout,*) ' '
    end if
  end subroutine domain_check

end module mod_clm_domain
