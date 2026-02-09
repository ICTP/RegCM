module mod_clm_subgrid
  !
  ! sub-grid data and mapping types and modules
  !
  use mod_realkinds
  use mod_intkinds
  use mod_clm_varpar
  use mod_clm_varctl
  use mod_clm_varsur

  implicit none (type, external)

  private

  save

  public :: subgrid_get_gcellinfo  ! Returns g,l,c,p properties from wtxy

  contains
  !
  ! Obtain gridcell properties
  !
  subroutine subgrid_get_gcellinfo (nw,                   &
                             nlunits, ncols, npfts,       &
                             nveg, wtveg,                 &
                             ncrop, wtcrop,               &
                             iurb, nurban_class, wturban_class, &
                             nlake, wtlake,               &
                             nwetland, wtwetland,         &
                             nglacier, wtglacier)
    implicit none (type, external)
    ! wtxy cell index
    integer(ik4), intent(in)  :: nw
    ! number of landunits
    integer(ik4), optional, intent(out) :: nlunits
    ! number of columns
    integer(ik4), optional, intent(out) :: ncols
    ! number of pfts
    integer(ik4), optional, intent(out) :: npfts
    ! number of vegetated pfts in naturally vegetated landunit
    integer(ik4), optional, intent(out) :: nveg
    ! weight (relative to gridcell) of naturally vegetated landunit
    real(rk8), optional, intent(out) :: wtveg
    ! number of crop pfts in crop landunit
    integer(ik4), optional, intent(out) :: ncrop
    ! weight (relative to gridcell) of crop landunit
    real(rk8), optional, intent(out) :: wtcrop
    ! Urban class index
    integer(ik4), optional, intent(in) :: iurb
    ! number of urban pfts (columns) in urban landunit
    integer(ik4), optional, intent(out) :: nurban_class
    ! weight (relative to gridcell) of urban pfts (columns) in urban landunit
    real(rk8), optional, intent(out) :: wturban_class
    ! number of lake pfts (columns) in lake landunit
    integer(ik4), optional, intent(out) :: nlake
    ! weight (relative to gridcell) of lake landunitof lake pfts (columns)
    ! in lake landunit
    real(rk8), optional, intent(out) :: wtlake
    ! number of wetland pfts (columns) in wetland landunit
    integer(ik4), optional, intent(out) :: nwetland
    ! weight (relative to gridcell) of wetland landunitof wetland pfts
    ! (columns) in wetland landunit
    real(rk8), optional, intent(out) :: wtwetland
    ! number of glacier pfts (columns) in glacier landunit
    integer(ik4), optional, intent(out) :: nglacier
    real(rk8), optional, intent(out) :: wtglacier
    integer(ik4) :: m, i             ! loop index
    integer(ik4) :: ipfts            ! number of pfts in gridcell
    integer(ik4) :: icols            ! number of columns in gridcell
    integer(ik4) :: ilunits          ! number of landunits in gridcell
    integer(ik4) :: npfts_per_lunit  ! number of pfts in landunit
    real(rk8) :: wtlunit             ! weight (relative to gridcell) of landunit

    ! Initialize pfts, columns and landunits counters for gridcell

    ipfts   = 0
    icols   = 0
    ilunits = 0

    ! Set naturally vegetated landunit

    npfts_per_lunit = 0
    wtlunit = 0._rk8
    ! If crop should be on separate land units
    if ( allocate_all_vegpfts .and. create_crop_landunit ) then
      do m = 1, maxpatch_pft-numcft
        if ( wtxy(nw,m) > 0.0_rk8 ) then
          npfts_per_lunit = npfts_per_lunit + 1 ! sum natural pfts
          wtlunit = wtlunit + wtxy(nw,m)        ! and their wts
        end if
      end do
      do m = maxpatch_pft-numcft+1, maxpatch_pft
        if ( wtxy(nw,m) > 0.0_rk8 ) then
          npfts_per_lunit = npfts_per_lunit + 1 ! sum crops, too, but not
        end if                                  ! their wts for now
      end do
    ! Assume that the vegetated landunit has one column
    else
      do m = 1, maxpatch_pft
        if ( wtxy(nw,m) > 0.0_rk8 ) then
          npfts_per_lunit = npfts_per_lunit + 1
          wtlunit = wtlunit + wtxy(nw,m)
        end if
      end do
    end if
    if ( npfts_per_lunit > 0 ) then ! true even when only crops are present
      if ( allocate_all_vegpfts ) then
        npfts_per_lunit = numpft+1
      end if
      if ( allocate_all_vegpfts .and. create_crop_landunit ) then
        npfts_per_lunit = numpft+1-numcft
      end if
      ilunits = ilunits + 1
      icols = icols + 1
    end if
    ipfts = ipfts + npfts_per_lunit
    if ( present(nveg ) ) nveg  = npfts_per_lunit
    if ( present(wtveg) ) wtveg = wtlunit

    ! Set urban landunit

    do i = 1, numurbl-1
      npfts_per_lunit = 0
      wtlunit = 0._rk8
      do m = npatch_urban_class(i), npatch_urban_class(i+1)-1
        if ( wtxy(nw,m) > 0.0_rk8 ) then
          npfts_per_lunit = npfts_per_lunit + 1
          wtlunit = wtlunit + wtxy(nw,m)
        end if
      end do
      if ( npfts_per_lunit > 0 ) then
        ilunits = ilunits + 1
        icols   = icols + npfts_per_lunit
      end if
      ipfts = ipfts + npfts_per_lunit
      if ( present(iurb) ) then
        if ( i == iurb ) then
          if ( present(nurban_class ) ) then
            nurban_class  = npfts_per_lunit
          end if
          if ( present(wturban_class) ) then
            wturban_class = wtlunit
          end if
        end if
      end if
    end do

    npfts_per_lunit = 0
    wtlunit = 0._rk8
    do m = npatch_urban_class(numurbl), npatch_lake-1
      if ( wtxy(nw,m) > 0.0_rk8 ) then
        npfts_per_lunit = npfts_per_lunit + 1
        wtlunit = wtlunit + wtxy(nw,m)
      end if
    end do
    if ( npfts_per_lunit > 0 ) then
      ilunits = ilunits + 1
      icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if ( present(iurb) ) then
      if ( iurb == numurbl ) then
        if ( present(nurban_class) ) then
          nurban_class = npfts_per_lunit
        end if
        if ( present(wturban_class) ) then
          wturban_class = wtlunit
        end if
      end if
    end if

    ! Set lake landunit

    npfts_per_lunit = 0
    wtlunit = 0._rk8
    if ( wtxy(nw,npatch_lake) > 0.0_rk8 ) then
      npfts_per_lunit = npfts_per_lunit + 1
      wtlunit = wtlunit + wtxy(nw,npatch_lake)
    end if
    if ( npfts_per_lunit > 0 ) then
      ilunits = ilunits + 1
      icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if ( present(nlake ) ) nlake  = npfts_per_lunit
    if ( present(wtlake) ) wtlake = wtlunit

    ! Set wetland landunit

    npfts_per_lunit = 0
    wtlunit = 0._rk8
    if ( wtxy(nw,npatch_wet) > 0.0_rk8 ) then
      npfts_per_lunit = npfts_per_lunit + 1
      wtlunit = wtlunit + wtxy(nw,npatch_wet)
    end if
    if ( npfts_per_lunit > 0 ) then
      ilunits = ilunits + 1
      icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if ( present(nwetland ) ) nwetland  = npfts_per_lunit
    if ( present(wtwetland) ) wtwetland = wtlunit

    ! Set glacier landunit

    npfts_per_lunit = 0
    wtlunit = 0._rk8
    if ( wtxy(nw,npatch_glacier) > 0.0_rk8 ) then
      npfts_per_lunit = npfts_per_lunit + 1
      wtlunit = wtlunit + wtxy(nw,npatch_glacier)
    end if
    if ( npfts_per_lunit > 0 ) then
      ilunits = ilunits + 1
      icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if ( present(nglacier ) ) nglacier  = npfts_per_lunit
    if ( present(wtglacier) ) wtglacier = wtlunit

    ! Set crop landunit if appropriate

    npfts_per_lunit = 0
    wtlunit = 0._rk8
    if ( allocate_all_vegpfts .and. create_crop_landunit ) then
      do m = 1, maxpatch_pft-numcft
        if ( wtxy(nw,m) > 0.0_rk8 ) then
          npfts_per_lunit = npfts_per_lunit + 1 ! sum natural pfts again
        end if                                   ! not their wts this time
      end do
      do m = maxpatch_pft-numcft+1, maxpatch_pft
        if ( wtxy(nw,m) > 0.0_rk8 ) then
          npfts_per_lunit = npfts_per_lunit + 1 ! sum crops
          wtlunit = wtlunit + wtxy(nw,m)        ! and their wts
        end if
      end do
    end if
    if ( npfts_per_lunit > 0 ) then ! true even if only natural veg is present
      if ( allocate_all_vegpfts .and. create_crop_landunit ) then
        npfts_per_lunit = numcft
      end if
      ilunits = ilunits + 1
      icols   = icols + npfts_per_lunit
    end if
    ipfts = ipfts + npfts_per_lunit
    if ( present(ncrop ) ) ncrop  = npfts_per_lunit
    if ( present(wtcrop) ) wtcrop = wtlunit

    ! Determine return arguments

    if ( present(nlunits) ) nlunits = ilunits
    if ( present(ncols) )   ncols   = icols
    if ( present(npfts) )   npfts   = ipfts

  end subroutine subgrid_get_gcellinfo

end module mod_clm_subgrid
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
