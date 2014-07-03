module mod_clm_cnverticalprofile
#ifdef CN
  !
  ! Module holding routines for vertical discretization of C and N
  ! inputs into deocmposing pools
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_mpmessage
  use mod_clm_varcon, only: dzsoi_decomp

  implicit none

  private

  save

  public:: decomp_vertprofiles

#ifdef VERTSOILC
  logical , public :: exponential_rooting_profile = .true.
  logical , public :: pftspecific_rootingprofile = .true.
  ! parameter for how steep the profile is for root C inputs
  ! (1/ e-folding depth) (1/m)
  real(rk8) , public :: rootprof_exp  = 3.0D0
  ! parameter for how steep the profile is for surface components
  ! (1/ e_folding depth) (1/m)
  real(rk8) , public :: surfprof_exp  = 10.0D0
#endif

  contains

  subroutine decomp_vertprofiles(lbp,ubp,lbc,ubc, &
                  num_soilc,filter_soilc,num_soilp,filter_soilp)
    use mod_clm_type
    use mod_clm_subgridave , only : p2c
    use mod_clm_varcon , only : zsoi, dzsoi, zisoi
    use mod_clm_varpar , only : nlevdecomp, nlevgrnd, &
            nlevdecomp_full, maxpatch_pft
    use mod_clm_pftvarcon , only : rootprof_beta, noveg
    implicit none
    integer(ik4), intent(in) :: lbp, ubp        ! pft-index bounds
    integer(ik4), intent(in) :: lbc, ubc        ! column-index bounds
    integer(ik4), intent(in) :: num_soilc       ! number of soil columns in filter
    integer(ik4), intent(in) :: filter_soilc(:) ! filter for soil columns
    integer(ik4), intent(in) :: num_soilp       ! number of soil pfts in filter
    integer(ik4), intent(in) :: filter_soilp(:) ! filter for soil pfts
    ! column level
    ! (1/m) profile for N fixation additions
    real(rk8), pointer :: nfixation_prof(:,:)
    ! (1/m) profile for N fixation additions
    real(rk8), pointer :: ndep_prof(:,:)
    integer(ik4), pointer :: altmax_lastyear_indx(:)  ! frost table depth (m)
    integer(ik4) , pointer :: npfts(:)  ! number of pfts for each column
    integer(ik4) , pointer :: pfti(:)   ! beginning pft index for each column

    ! pft level
    integer(ik4) , pointer :: ivt(:)    ! pft vegetation type
    ! fraction of roots in each soil layer  (nlevgrnd)
    real(rk8), pointer :: rootfr(:,:)
    integer(ik4) , pointer :: pcolumn(:)    ! pft's column index
    real(rk8), pointer :: leaf_prof(:,:)    ! (1/m) profile of leaves
    real(rk8), pointer :: froot_prof(:,:)   ! (1/m) profile of fine roots
    real(rk8), pointer :: croot_prof(:,:)   ! (1/m) profile of coarse roots
    real(rk8), pointer :: stem_prof(:,:)    ! (1/m) profile of stems
    real(rk8), pointer :: wtcol(:)          ! pft weight relative to column (0-1)
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer :: pactive(:)

    ! local variables
    real(rk8) :: surface_prof(1:nlevdecomp)
    ! pft-native root fraction used for calculating inputs
    real(rk8) :: cinput_rootfr(lbp:ubp, 1:nlevdecomp_full)
#ifdef VERTSOILC
    real(rk8) :: surface_prof_tot
    real(rk8) :: rootfr_tot
    ! col-native root fraction used for calculating inputs
    real(rk8) :: col_cinput_rootfr(lbc:ubc, 1:nlevdecomp_full)
    integer(ik4) :: pi
#endif
    integer(ik4)  :: c, j, fc, p, fp

    ! debugging temp variables
    real(rk8) :: froot_prof_sum
    real(rk8) :: croot_prof_sum
    real(rk8) :: leaf_prof_sum
    real(rk8) :: stem_prof_sum
    real(rk8) :: ndep_prof_sum
    real(rk8) :: nfixation_prof_sum
    real(rk8) :: delta = 1.e-10

    ! assign local pointers at the column level
    nfixation_prof         => clm3%g%l%c%cps%nfixation_prof
    ndep_prof              => clm3%g%l%c%cps%ndep_prof
    altmax_lastyear_indx   => clm3%g%l%c%cps%altmax_lastyear_indx
    npfts                  => clm3%g%l%c%npfts
    pfti                   => clm3%g%l%c%pfti

    ! assign local pointers at the pft level
    ivt                    => clm3%g%l%c%p%itype
    leaf_prof              => clm3%g%l%c%p%pps%leaf_prof
    froot_prof             => clm3%g%l%c%p%pps%froot_prof
    croot_prof             => clm3%g%l%c%p%pps%croot_prof
    stem_prof              => clm3%g%l%c%p%pps%stem_prof
    pcolumn                => clm3%g%l%c%p%column
    rootfr                 => clm3%g%l%c%p%pps%rootfr
    wtcol                  => clm3%g%l%c%p%wtcol
    pactive                => clm3%g%l%c%p%active


#ifdef VERTSOILC
    ! define a single shallow surface profile for surface additions
    ! (leaves, stems, and N deposition)
    surface_prof(:) = 0.D0
    do j = 1 , nlevdecomp
      surface_prof(j) = exp(-surfprof_exp * zsoi(j)) / dzsoi_decomp(j)
    end do

    ! initialize profiles to zero
    leaf_prof(:,:) = 0.D0
    froot_prof(:,:) = 0.D0
    croot_prof(:,:) = 0.D0
    stem_prof(:,:) = 0.D0
    nfixation_prof(:,:) = 0.D0
    ndep_prof(:,:) = 0.D0

    cinput_rootfr(:,:) = 0.D0
    col_cinput_rootfr(:,:) = 0.D0

    if ( exponential_rooting_profile ) then
      if ( .not. pftspecific_rootingprofile ) then
        ! define rooting profile from exponential parameters
        do j = 1, nlevdecomp
          do fp = 1,num_soilp
            p = filter_soilp(fp)
            cinput_rootfr(p,j) = exp(-rootprof_exp * zsoi(j)) / dzsoi_decomp(j)
          end do
        end do
      else
        ! use beta distribution parameter from Jackson et al., 1996
        do p = lbp, ubp
          c = pcolumn(p)
          if (ivt(p) /= noveg) then
            do j = 1, nlevdecomp
              cinput_rootfr(p,j) = ( rootprof_beta(ivt(p)) ** &
                (zisoi(j-1)*100.D0) - rootprof_beta(ivt(p)) ** &
                (zisoi(j)*100.D0) ) / dzsoi_decomp(j)
            end do
          else
            cinput_rootfr(p,1) = 1.D0 / dzsoi_decomp(1)
          end if
        end do
      end if
    else
      do j = 1, nlevdecomp
        ! use standard CLM root fraction profiles
        do fp = 1,num_soilp
          p = filter_soilp(fp)
          cinput_rootfr(p,j) = rootfr(p,j) / dzsoi_decomp(j)
        end do
      end do
    end if

    do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)
      ! integrate rootfr over active layer of soil column
      rootfr_tot = 0.D0
      surface_prof_tot = 0.D0
      do j = 1, min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
        rootfr_tot = rootfr_tot + cinput_rootfr(p,j) * dzsoi_decomp(j)
        surface_prof_tot = surface_prof_tot + surface_prof(j)  * dzsoi_decomp(j)
      end do
      if ( (altmax_lastyear_indx(c) .gt. 0) .and. &
           (rootfr_tot .gt. 0.D0) .and. (surface_prof_tot .gt. 0.D0) ) then
        ! where there is not permafrost extending to the surface, integrate
        ! the profiles over the active layer
        ! this is equivalent to integrating over all soil layers outside
        ! of permafrost regions
        do j = 1, min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
          froot_prof(p,j) = cinput_rootfr(p,j) / rootfr_tot
          croot_prof(p,j) = cinput_rootfr(p,j) / rootfr_tot
          ! set all surface processes to shallower profile
          leaf_prof(p,j) = surface_prof(j)/ surface_prof_tot
          stem_prof(p,j) = surface_prof(j)/ surface_prof_tot
        end do
      else
        ! if fully frozen, or no roots, put everything in the top layer
        froot_prof(p,1) = 1.0D0/dzsoi_decomp(1)
        croot_prof(p,1) = 1.0D0/dzsoi_decomp(1)
        leaf_prof(p,1) = 1.0D0/dzsoi_decomp(1)
        stem_prof(p,1) = 1.0D0/dzsoi_decomp(1)
      end if
    end do

    !! aggregate root profile to column
    ! call p2c (lbp, ubp, lbc, ubc, nlevdecomp_full, cinput_rootfr,
    !           col_cinput_rootfr, 'unity')
    do pi = 1,maxpatch_pft
      do fc = 1,num_soilc
        c = filter_soilc(fc)
        if (pi <=  npfts(c)) then
          p = pfti(c) + pi - 1
          if (pactive(p)) then
            do j = 1,nlevdecomp
              col_cinput_rootfr(c,j) = col_cinput_rootfr(c,j) + &
                      cinput_rootfr(p,j) * wtcol(p)
            end do
          end if
        end if
      end do
    end do

    ! repeat for column-native profiles: Ndep and Nfix
    do fc = 1,num_soilc
      c = filter_soilc(fc)
      rootfr_tot = 0.D0
      surface_prof_tot = 0.D0
      ! redo column ntegration over active layer for column-native profiles
      do j = 1, min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
        rootfr_tot = rootfr_tot + col_cinput_rootfr(c,j) * dzsoi_decomp(j)
        surface_prof_tot = surface_prof_tot + surface_prof(j) * dzsoi_decomp(j)
      end do
      if ( (altmax_lastyear_indx(c) .gt. 0) .and. &
           (rootfr_tot .gt. 0.D0) .and. (surface_prof_tot .gt. 0.D0) ) then
        do j = 1,  min(max(altmax_lastyear_indx(c), 1), nlevdecomp)
          nfixation_prof(c,j) = col_cinput_rootfr(c,j) / rootfr_tot
          ndep_prof(c,j) = surface_prof(j)/ surface_prof_tot
        end do
      else
        nfixation_prof(c,1) = 1./dzsoi_decomp(1)
        ndep_prof(c,1) = 1./dzsoi_decomp(1)
      end if
    end do

#else

    ! for one layer decomposition model, set profiles to unity
    leaf_prof(:,:) = 1.D0
    froot_prof(:,:) = 1.D0
    croot_prof(:,:) = 1.D0
    stem_prof(:,:) = 1.D0
    nfixation_prof(:,:) = 1.D0
    ndep_prof(:,:) = 1.D0

#endif


    ! check to make sure integral of all profiles = 1.
    do fc = 1,num_soilc
      c = filter_soilc(fc)
      ndep_prof_sum = 0.
      nfixation_prof_sum = 0.
      do j = 1, nlevdecomp
        ndep_prof_sum = ndep_prof_sum + ndep_prof(c,j) *  dzsoi_decomp(j)
        nfixation_prof_sum = nfixation_prof_sum + &
                nfixation_prof(c,j) *  dzsoi_decomp(j)
      end do
      if ( ( abs(ndep_prof_sum - 1.D0) .gt. delta ) .or. &
           ( abs(nfixation_prof_sum - 1.D0) .gt. delta ) ) then
        write(stderr, *) 'profile sums: ', ndep_prof_sum, nfixation_prof_sum
        write(stderr, *) 'c: ', c
        write(stderr, *) 'altmax_lastyear_indx: ', altmax_lastyear_indx(c)
        write(stderr, *) 'nfixation_prof: ', nfixation_prof(c,:)
        write(stderr, *) 'ndep_prof: ', ndep_prof(c,:)
        write(stderr, *) 'cinput_rootfr: ', cinput_rootfr(c,:)
        write(stderr, *) 'dzsoi_decomp: ', dzsoi_decomp(:)
        write(stderr, *) 'surface_prof: ', surface_prof(:)
        write(stderr, *) 'npfts(c): ', npfts(c)
        do p = pfti(c), pfti(c) + npfts(c) -1
          write(stderr, *) 'p, ivt(p), wtcol(p): ', p, ivt(p), wtcol(p)
          write(stderr, *) 'cinput_rootfr(p,:): ', cinput_rootfr(p,:)
        end do
        call fatal(__FILE__,__LINE__, &
            "decomp_vertprofiles ERROR: _prof_sum-1>delta" )
      end if
    end do

    do fp = 1,num_soilp
      p = filter_soilp(fp)
      froot_prof_sum = 0.0D0
      croot_prof_sum = 0.0D0
      leaf_prof_sum = 0.0D0
      stem_prof_sum = 0.0D0
      do j = 1, nlevdecomp
        froot_prof_sum = froot_prof_sum + froot_prof(p,j) *  dzsoi_decomp(j)
        croot_prof_sum = croot_prof_sum + croot_prof(p,j) *  dzsoi_decomp(j)
        leaf_prof_sum = leaf_prof_sum + leaf_prof(p,j) *  dzsoi_decomp(j)
        stem_prof_sum = stem_prof_sum + stem_prof(p,j) *  dzsoi_decomp(j)
      end do
      if ( ( abs(froot_prof_sum - 1.D0) .gt. delta ) .or. &
           ( abs(croot_prof_sum - 1.D0) .gt. delta ) .or. &
           ( abs(stem_prof_sum - 1.D0) .gt. delta ) .or.  &
           ( abs(leaf_prof_sum - 1.D0) .gt. delta ) ) then
        write(stderr, *) 'profile sums: ', froot_prof_sum, &
                croot_prof_sum, leaf_prof_sum, stem_prof_sum
        call fatal(__FILE__,__LINE__, &
            'decomp_vertprofiles ERROR: sum-1 > delta' )
      end if
    end do
  end subroutine decomp_vertprofiles

#endif

end module mod_clm_cnverticalprofile
