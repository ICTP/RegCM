module mod_clm_cndvecosystemdynini

#if (defined CNDV)
  use mod_intkinds
  use mod_realkinds

  implicit none

  private

  save
  public  :: CNDVEcosystemDynini ! CNDV related initializations

  contains
  !
  ! CNDV related initializations
  !
  subroutine CNDVEcosystemDynini(adomain)
    use mod_clm_type
    use mod_clm_decomp , only : get_proc_bounds , get_proc_global
    use mod_constants , only : tzero
    use mod_clm_atmlnd , only : atm_domain
    implicit none
    type(atm_domain) , intent(in) :: adomain
    integer(ik4) :: g , p , n   ! indices
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc !              "                column indices
    integer(ik4) :: begl , endl !              "                landunit indices
    integer(ik4) :: begg , endg !              "                gridcell indices
    type(gridcell_type) , pointer :: gptr ! pointer to gridcell derived subtype
    type(pft_type) , pointer :: pptr      ! pointer to pft derived subtype

    ! Set pointers into derived type

    gptr => clm3%g
    pptr => clm3%g%l%c%p

    ! ---------------------------------------------------------------
    ! Some of the following came from LPJ subroutine initgrid
    ! ---------------------------------------------------------------

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    do p = begp , endp
      pptr%pdgvs%present(p)   = .false.
      pptr%pdgvs%crownarea(p) = 0._rkx
      pptr%pdgvs%nind(p)      = 0._rkx
      pptr%pcs%leafcmax(p)    = 0._rkx
      pptr%pdgvs%t_mo_min(p)  = 1.0e+36_rkx
    end do

    do g = begg , endg
      gptr%gdgvs%agdd20(g)   = 0._rkx
      !initialize this way for Phenology code
      gptr%gdgvs%tmomin20(g) = adomain%tgrd(g) - 5._rkx
    end do
  end subroutine CNDVEcosystemDynini

#endif

end module mod_clm_cndvecosystemdynini
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
