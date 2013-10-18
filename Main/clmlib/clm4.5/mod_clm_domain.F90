module mod_clm_domain
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: domainMod
!
! !DESCRIPTION:
! Module containing 2-d global surface boundary data information
!
! !USES:
  use mod_realkinds
  use mod_mppparam
  use mod_runparams
  use mod_mpmessage
!
! !PUBLIC TYPES:
  implicit none
  private
!
  public :: domain_type

  !--- this typically contains local domain info with arrays dim begg:endg ---
  type domain_type
     integer          :: ns         ! global size of domain
     integer          :: ni,nj      ! global axis if 2d (nj=1 if unstructured)
     logical          :: isgrid2d   ! true => global grid is lat/lon
     integer          :: nbeg,nend  ! local beg/end indices
     character(len=8) :: clmlevel   ! grid type
     integer ,pointer :: mask(:)    ! land mask: 1 = land, 0 = ocean
     real(rk8),pointer :: frac(:)    ! fractional land
     real(rk8),pointer :: topo(:)    ! topography
     real(rk8),pointer :: latc(:)    ! latitude of grid cell (deg)
     real(rk8),pointer :: lonc(:)    ! longitude of grid cell (deg)
     real(rk8),pointer :: area(:)    ! grid cell area (km**2)
     integer ,pointer :: pftm(:)    ! pft mask: 1=real, 0=fake, -1=notset
     integer ,pointer :: glcmask(:) ! glc mask: 1=sfc mass balance required by GLC component
                                    ! 0=SMB not required (default)
     character(len=16):: set        ! flag to check if domain is set
     logical          :: decomped   ! decomposed locally or global copy
  end type domain_type

  type(domain_type)    , public :: ldomain
  real(rk8), allocatable, public :: lon1d(:), lat1d(:) ! 1d lat/lons for 2d grids
!
! !PUBLIC MEMBER FUNCTIONS:
  public domain_init          ! allocates/nans domain types
  public domain_clean         ! deallocates domain types
  public domain_check         ! write out domain info
!
! !REVISION HISTORY:
! Originally clm_varsur by Mariana Vertenstein
! Migrated from clm_varsur to domainMod by T Craig
!
  character(len=16),parameter :: set   = 'domain_set      '
  character(len=16),parameter :: unset = 'NOdomain_unsetNO'
!
!EOP
!------------------------------------------------------------------------------

contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_init
!
! !INTERFACE:
  subroutine domain_init(domain,isgrid2d,ni,nj,nbeg,nend,clmlevel)
!
! !DESCRIPTION:
! This subroutine allocates and nans the domain type
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type)   :: domain        ! domain datatype
    logical, intent(in) :: isgrid2d      ! true => global grid is lat/lon
    integer, intent(in) :: ni,nj         ! grid size, 2d
    integer         , intent(in), optional  :: nbeg,nend  ! beg/end indices
    character(len=*), intent(in), optional  :: clmlevel   ! grid type
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer ier
    integer nb,ne
!
!------------------------------------------------------------------------------

    nb = 1
    ne = ni*nj
    if (present(nbeg)) then
       if (present(nend)) then
          nb = nbeg
          ne = nend
       endif
    endif

    if (domain%set == set) then
       call domain_clean(domain)
    endif
    allocate(domain%mask(nb:ne),domain%frac(nb:ne),domain%latc(nb:ne), &
             domain%pftm(nb:ne),domain%area(nb:ne),domain%lonc(nb:ne), &
             domain%topo(nb:ne),domain%glcmask(nb:ne),stat=ier)
    if (ier /= 0) then
       call fatal(__FILE__,__LINE__, &
         'domain_init ERROR: allocate mask, frac, lat, lon, area ')
    endif

    if (present(clmlevel)) then
       domain%clmlevel = clmlevel
    endif

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
    if (domain%nbeg == 1 .and. domain%nend == domain%ns) then
       domain%decomped = .false.
    else
       domain%decomped = .true.
    endif

    domain%pftm     = -9999
    domain%glcmask  = 0  

end subroutine domain_init
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_clean
!
! !INTERFACE:
  subroutine domain_clean(domain)
!
! !DESCRIPTION:
! This subroutine deallocates the domain type
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type) :: domain        ! domain datatype
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!EOP
    integer ier
!
!------------------------------------------------------------------------------
    if (domain%set == set) then
       if (myid == italk) then
          write(stdout,*) 'domain_clean: cleaning ',domain%ni,domain%nj
       endif
       deallocate(domain%mask,domain%frac,domain%latc, &
                  domain%lonc,domain%area,domain%pftm, &
                  domain%topo,domain%glcmask,stat=ier)
       if (ier /= 0) then
          call fatal(__FILE__,__LINE__, &
            'domain_clean ERROR: deallocate mask, frac, lat, lon, area ')
       endif
    else
       if (myid == italk) then
          write(stdout,*) 'domain_clean WARN: clean domain unecessary '
       endif
    endif

    domain%clmlevel   = unset
    domain%ns         = bigint
    domain%ni         = bigint
    domain%nj         = bigint
    domain%nbeg       = bigint
    domain%nend       = bigint
    domain%set        = unset
    domain%decomped   = .true.

end subroutine domain_clean
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: domain_check
!
! !INTERFACE:
  subroutine domain_check(domain)
!
! !DESCRIPTION:
! This subroutine write domain info
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    type(domain_type),intent(in)  :: domain        ! domain datatype
!
! !REVISION HISTORY:
!   Created by T Craig
!
!
! !LOCAL VARIABLES:
!
!EOP
!------------------------------------------------------------------------------

  if (myid == italk) then
    write(stdout,*) '  domain_check set       = ',trim(domain%set)
    write(stdout,*) '  domain_check decomped  = ',domain%decomped
    write(stdout,*) '  domain_check ns        = ',domain%ns
    write(stdout,*) '  domain_check ni,nj     = ',domain%ni,domain%nj
    write(stdout,*) '  domain_check clmlevel  = ',trim(domain%clmlevel)
    write(stdout,*) '  domain_check nbeg,nend = ',domain%nbeg,domain%nend
    write(stdout,*) '  domain_check lonc      = ',minval(domain%lonc),maxval(domain%lonc)
    write(stdout,*) '  domain_check latc      = ',minval(domain%latc),maxval(domain%latc)
    write(stdout,*) '  domain_check mask      = ',minval(domain%mask),maxval(domain%mask)
    write(stdout,*) '  domain_check frac      = ',minval(domain%frac),maxval(domain%frac)
    write(stdout,*) '  domain_check topo      = ',minval(domain%topo),maxval(domain%topo)
    write(stdout,*) '  domain_check area      = ',minval(domain%area),maxval(domain%area)
    write(stdout,*) '  domain_check pftm      = ',minval(domain%pftm),maxval(domain%pftm)
    write(stdout,*) '  domain_check glcmask   = ',minval(domain%glcmask),maxval(domain%glcmask)
    write(stdout,*) ' '
  endif

end subroutine domain_check

!------------------------------------------------------------------------------

end module mod_clm_domain
