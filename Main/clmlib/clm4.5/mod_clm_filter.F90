module mod_clm_filter
  !
  ! Module of filters used for processing columns and pfts of particular
  ! types, including lake, non-lake, urban, soil, snow, non-snow, and
  ! naturally-vegetated patches.
  !
  use mod_realkinds
  use mod_intkinds
  use mod_mpmessage
  use mod_stdio
  use mod_clm_type
  use mod_clm_decomp , only : get_proc_bounds
  use mod_clm_pftvarcon , only : npcropmin
  use mod_clm_varcon , only : istsoil , isturb , icol_road_perv , istcrop
  use mod_dynparam

  implicit none

  private

  save

  type procfilter
#if (defined CNDV)
    ! CNDV nat-vegetated (present) filter (pfts)
    integer(ik4) , pointer , dimension(:) :: natvegp
    ! number of pfts in nat-vegetated filter
    integer(ik4) :: num_natvegp
#endif
    ! prognostic crop filter (pfts)
    integer(ik4) , pointer , dimension(:) :: pcropp
    ! number of pfts in prognostic crop filter
    integer(ik4) :: num_pcropp
    ! soil w/o prog. crops (pfts)
    integer(ik4) , pointer , dimension(:) :: soilnopcropp
    ! number of pfts in soil w/o prog crops
    integer(ik4) :: num_soilnopcropp

    ! lake filter (pfts)
    integer(ik4) , pointer , dimension(:) :: lakep
    ! number of pfts in lake filter
    integer(ik4) :: num_lakep
    ! non-lake filter (pfts)
    integer(ik4) , pointer , dimension(:) :: nolakep
    ! number of pfts in non-lake filter
    integer(ik4) :: num_nolakep
    ! lake filter (columns)
    integer(ik4) , pointer , dimension(:) :: lakec
    ! number of columns in lake filter
    integer(ik4) :: num_lakec
    ! non-lake filter (columns)
    integer(ik4) , pointer , dimension(:) :: nolakec
    ! number of columns in non-lake filter
    integer(ik4) :: num_nolakec

    ! soil filter (columns)
    integer(ik4) , pointer , dimension(:) :: soilc
    ! number of columns in soil filter
    integer(ik4) :: num_soilc
    ! soil filter (pfts)
    integer(ik4) , pointer , dimension(:) :: soilp
    ! number of pfts in soil filter
    integer(ik4) :: num_soilp

    ! snow filter (columns)
    integer(ik4) , pointer , dimension(:) :: snowc
    ! number of columns in snow filter
    integer(ik4) :: num_snowc
    ! non-snow filter (columns)
    integer(ik4) , pointer , dimension(:) :: nosnowc
    ! number of columns in non-snow filter
    integer(ik4) :: num_nosnowc

    ! hydrology filter (columns)
    integer(ik4) , pointer , dimension(:) :: hydrologyc
    ! number of columns in hydrology filter
    integer(ik4) :: num_hydrologyc

    ! urban filter (landunits)
    integer(ik4) , pointer , dimension(:) :: urbanl
    ! number of landunits in urban filter
    integer(ik4) :: num_urbanl
    ! non-urban filter (landunits)
    integer(ik4) , pointer , dimension(:) :: nourbanl
    ! number of landunits in non-urban filter
    integer(ik4) :: num_nourbanl

    ! urban filter (columns)
    integer(ik4) , pointer , dimension(:) :: urbanc
    ! number of columns in urban filter
    integer(ik4) :: num_urbanc
    ! non-urban filter (columns)
    integer(ik4) , pointer , dimension(:) :: nourbanc
    ! number of columns in non-urban filter
    integer(ik4) :: num_nourbanc

    ! urban filter (pfts)
    integer(ik4) , pointer , dimension(:) :: urbanp
    ! number of pfts in urban filter
    integer(ik4) :: num_urbanp
    ! non-urban filter (pfts)
    integer(ik4) , pointer , dimension(:) :: nourbanp
    ! number of pfts in non-urban filter
    integer(ik4) :: num_nourbanp

    ! non-lake, non-urban filter (pfts)
    integer(ik4) , pointer , dimension(:) :: nolakeurbanp
    ! number of pfts in non-lake, non-urban filter
    integer(ik4) :: num_nolakeurbanp
  end type procfilter

  public :: procfilter

  type(procfilter), public :: filter

  public :: allocFilters   ! allocate memory for filters
  public :: setFilters     ! set filters

  contains
  !
  ! Allocate CLM filters.
  !
  subroutine allocFilters()
    implicit none
    integer(ik4) :: begp , endp  ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc  ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl  ! per-proc beginning and ending ldunit indices
    integer(ik4) :: begg , endg  ! per-proc beginning and ending gdcell indices
    integer(ik4) :: ier          ! error status

    ! Determine variables for this processor

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    allocate(filter%lakep(endp-begp+1))
    allocate(filter%nolakep(endp-begp+1))
    allocate(filter%nolakeurbanp(endp-begp+1))

    allocate(filter%lakec(endc-begc+1))
    allocate(filter%nolakec(endc-begc+1))

    allocate(filter%soilc(endc-begc+1))
    allocate(filter%soilp(endp-begp+1))

    allocate(filter%snowc(endc-begc+1))
    allocate(filter%nosnowc(endc-begc+1))

#if (defined CNDV)
    allocate(filter%natvegp(endp-begp+1))
#endif

    allocate(filter%hydrologyc(endc-begc+1))

    allocate(filter%urbanp(endp-begp+1))
    allocate(filter%nourbanp(endp-begp+1))

    allocate(filter%urbanc(endc-begc+1))
    allocate(filter%nourbanc(endc-begc+1))

    allocate(filter%urbanl(endl-begl+1))
    allocate(filter%nourbanl(endl-begl+1))

    allocate(filter%pcropp(endp-begp+1))
    allocate(filter%soilnopcropp(endp-begp+1))

    filter%lakep = -1
    filter%nolakep = -1
    filter%nolakeurbanp = -1
    filter%lakec = -1
    filter%nolakec = -1
    filter%soilc = -1
    filter%soilp = -1
    filter%snowc = -1
    filter%nosnowc = -1
#if (defined CNDV)
    filter%natvegp = -1
#endif
    filter%hydrologyc = -1
    filter%urbanp = -1
    filter%nourbanp = -1
    filter%urbanc = -1
    filter%nourbanc = -1
    filter%urbanl = -1
    filter%nourbanl = -1
    filter%pcropp = -1
    filter%soilnopcropp = -1
  end subroutine allocFilters
  !
  ! Set CLM filters.
  !
  subroutine setFilters( )
    implicit none
    integer(ik4) , pointer , dimension(:) :: ctype ! column type
    integer(ik4) :: c , l , p   ! column, landunit, pft indices
    integer(ik4) :: fl          ! lake filter index
    integer(ik4) :: fnl , fnlu  ! non-lake filter index
    integer(ik4) :: fs          ! soil filter index
    integer(ik4) :: f , fn      ! general indices
    integer(ik4) :: begp , endp  ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc  ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl  ! per-proc beginning and ending ldunit indices
    integer(ik4) :: begg , endg  ! per-proc beginning and ending gdcell indices

    ctype => clm3%g%l%c%itype

    ! Determine boundaries

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! Create lake and non-lake filters at column-level

    fl = 0
    fnl = 0
    do c = begc , endc
      l = clm3%g%l%c%landunit(c)
      if ( clm3%g%l%lakpoi(l) ) then
        fl = fl + 1
        filter%lakec(fl) = c
      else
        fnl = fnl + 1
        filter%nolakec(fnl) = c
      end if
    end do
    filter%num_lakec = fl
    filter%num_nolakec = fnl

    ! Create lake and non-lake filters at pft-level
    ! Filter will only be active if pft is active

    fl = 0
    fnl = 0
    fnlu = 0
    do p = begp , endp
      if ( clm3%g%l%c%p%active(p) ) then
        l = clm3%g%l%c%p%landunit(p)
        if ( clm3%g%l%lakpoi(l) ) then
          fl = fl + 1
          filter%lakep(fl) = p
        else
          fnl = fnl + 1
          filter%nolakep(fnl) = p
          if ( clm3%g%l%itype(l) /= isturb ) then
            fnlu = fnlu + 1
            filter%nolakeurbanp(fnlu) = p
          end if
        end if
      end if
    end do
    filter%num_lakep = fl
    filter%num_nolakep = fnl
    filter%num_nolakeurbanp = fnlu

    ! Create soil filter at column-level

    fs = 0
    do c = begc , endc
      l = clm3%g%l%c%landunit(c)
      if ( clm3%g%l%itype(l) == istsoil .or. &
           clm3%g%l%itype(l) == istcrop ) then
        fs = fs + 1
        filter%soilc(fs) = c
      end if
    end do
    filter%num_soilc = fs

    ! Create soil filter at pft-level
    ! Filter will only be active if pft is active

    fs = 0
    do p = begp , endp
      if ( clm3%g%l%c%p%active(p) ) then
        l = clm3%g%l%c%p%landunit(p)
        if ( clm3%g%l%itype(l) == istsoil .or. &
             clm3%g%l%itype(l) == istcrop ) then
          fs = fs + 1
          filter%soilp(fs) = p
        end if
      end if
    end do
    filter%num_soilp = fs

    ! Create column-level hydrology filter (soil and Urban pervious road cols)

    f = 0
    do c = begc , endc
      l = clm3%g%l%c%landunit(c)
      if ( clm3%g%l%itype(l) == istsoil .or. &
           ctype(c) == icol_road_perv .or.   &
           clm3%g%l%itype(l) == istcrop ) then
        f = f + 1
        filter%hydrologyc(f) = c
      end if
    end do
    filter%num_hydrologyc = f

    ! Create prognostic crop and soil w/o prog. crop filters at pft-level
    ! according to where the crop model should be used

    fl  = 0
    fnl = 0
    do p = begp , endp
      if ( clm3%g%l%c%p%active(p) ) then
        if ( clm3%g%l%c%p%itype(p) >= npcropmin ) then
          ! skips 2 generic crop types
          fl = fl + 1
          filter%pcropp(fl) = p
        else
          l = clm3%g%l%c%p%landunit(p)
          if ( clm3%g%l%itype(l) == istsoil .or. &
               clm3%g%l%itype(l) == istcrop ) then
            fnl = fnl + 1
            filter%soilnopcropp(fnl) = p
          end if
        end if
      end if
    end do
    filter%num_pcropp   = fl
    filter%num_soilnopcropp = fnl   ! This wasn't being set before...

    ! Create landunit-level urban and non-urban filters

    f = 0
    fn = 0
    do l = begl , endl
      if ( clm3%g%l%itype(l) == isturb ) then
        f = f + 1
        filter%urbanl(f) = l
      else
        fn = fn + 1
        filter%nourbanl(fn) = l
      end if
    end do
    filter%num_urbanl = f
    filter%num_nourbanl = fn

    ! Create column-level urban and non-urban filters

    f = 0
    fn = 0
    do c = begc , endc
      l = clm3%g%l%c%landunit(c)
      if ( clm3%g%l%itype(l) == isturb ) then
        f = f + 1
        filter%urbanc(f) = c
      else
        fn = fn + 1
        filter%nourbanc(fn) = c
      end if
    end do
    filter%num_urbanc = f
    filter%num_nourbanc = fn

    ! Create pft-level urban and non-urban filters

    f = 0
    fn = 0
    do p = begp , endp
      l = clm3%g%l%c%p%landunit(p)
      if ( clm3%g%l%itype(l) == isturb .and. clm3%g%l%c%p%active(p) ) then
        f = f + 1
        filter%urbanp(f) = p
      else
        fn = fn + 1
        filter%nourbanp(fn) = p
      end if
    end do
    filter%num_urbanp = f
    filter%num_nourbanp = fn

    ! Note: snow filters are reconstructed each time step in Hydrology2
    ! Note: CNDV "pft present" filter is reconstructed each time CNDV is run

  end subroutine setFilters

end module mod_clm_filter
