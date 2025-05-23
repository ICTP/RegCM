module mod_clm_staticecosysdyn

  ! Static Ecosystem dynamics: phenology, vegetation.
  ! This is for the CLM Satelitte Phenology model (CLMSP). Allow some
  ! subroutines to be used by the CLM Carbon Nitrogen model (CLMCN)
  ! so that DryDeposition code can get estimates of LAI differences
  ! between months.

  use mod_realkinds
  use mod_intkinds
  use mod_mpmessage
  use mod_stdio
  use mod_date
  use mod_dynparam
  use mod_mppparam
  use mod_runparams
  use mod_mpmessage
  use mod_clm_nchelper
  use mod_clm_decomp, only : get_proc_bounds, gcomm_gridcell
  use mod_clm_type
  use mod_clm_pftvarcon, only : noveg, nc3crop, nbrdlf_dcd_brl_shrub
  use mod_clm_varctl, only : fsurdat, nextdate
  use mod_clm_varpar, only : numpft
  use mod_clm_domain, only : ldomain

  implicit none

  private

  save

#if (!defined CN) && (!defined CNDV)
  ! CLMSP Ecosystem dynamics: phenology, vegetation
  public :: EcosystemDyn
#endif
  ! Dynamically allocate memory
  public :: EcosystemDynini
  ! interpolate monthly vegetation data
  public :: interpMonthlyVeg
  ! Read in annual vegetation (needed for Dry-deposition)
  public :: readAnnualVegetation

  ! read monthly vegetation data for two months
  private :: readMonthlyVegetation

  integer(ik4) :: InterpMonths1 = -9999 ! saved month index
  real(rk8), dimension(2) :: timwt  ! time weights for month 1 and month 2
  ! lai for interpolation (2 months)
  real(rk8), allocatable, dimension(:,:) :: mlai2t
  ! sai for interpolation (2 months)
  real(rk8), allocatable, dimension(:,:) :: msai2t
  ! top vegetation height for interpolation (2 months)
  real(rk8), allocatable, dimension(:,:) :: mhvt2t
  ! bottom vegetation height for interpolation(2 months)
  real(rk8), allocatable, dimension(:,:) :: mhvb2t

  contains
  !
  ! Dynamically allocate memory and set to signaling NaN.
  !
  subroutine EcosystemDynIni()
    implicit none
    integer(ik4) :: ier          ! error code
    integer(ik4) :: begp, endp  ! local beg and end p index

    call get_proc_bounds(begp=begp,endp=endp)
    ier = 0
    if(.not.allocated(mlai2t))allocate (mlai2t(begp:endp,2), &
              msai2t(begp:endp,2), &
              mhvt2t(begp:endp,2), &
              mhvb2t(begp:endp,2), stat=ier)
    if (ier /= 0) then
       write(stderr,*) 'EcosystemDynini allocation error'
       call fatal(__FILE__,__LINE__,'clm now stopping')
    end if
    mlai2t(:,:) = nan
    msai2t(:,:) = nan
    mhvt2t(:,:) = nan
    mhvb2t(:,:) = nan
  end subroutine EcosystemDynini

#if (!defined CN) && (!defined CNDV)
  !
  ! Ecosystem dynamics: phenology, vegetation
  ! Calculates leaf areas (tlai, elai),  stem areas (tsai, esai) and
  ! height (htop).
  !
  subroutine EcosystemDyn(lbp, ubp, num_nolakep, filter_nolakep, doalb)
    implicit none
    ! pft bounds
    integer(ik4), intent(in) :: lbp, ubp
    ! number of column non-lake points in pft filter
    integer(ik4), intent(in) :: num_nolakep
    ! pft filter for non-lake points
    integer(ik4), intent(in), dimension(ubp-lbp+1) :: filter_nolakep
    ! true = surface albedo calculation time step
    logical, intent(in) :: doalb
    ! fraction of ground covered by snow (0 to 1)
    real(rk8), pointer, contiguous, dimension(:) :: frac_sno
    ! column index associated with each pft
    integer(ik4), pointer, contiguous, dimension(:) :: pcolumn
    ! snow height (m)
    real(rk8), pointer, contiguous, dimension(:) :: snow_depth
    ! pft vegetation type
    integer(ik4), pointer, contiguous, dimension(:) :: ivt
    ! one-sided leaf area index, no burying by snow
    real(rk8), pointer, contiguous, dimension(:) :: tlai
    ! one-sided stem area index, no burying by snow
    real(rk8), pointer, contiguous, dimension(:) :: tsai
    ! canopy top (m)
    real(rk8), pointer, contiguous, dimension(:) :: htop
    ! canopy bottom (m)
    real(rk8), pointer, contiguous, dimension(:) :: hbot
    ! one-sided leaf area index with burying by snow
    real(rk8), pointer, contiguous, dimension(:) :: elai
    ! one-sided stem area index with burying by snow
    real(rk8), pointer, contiguous, dimension(:) :: esai
    ! frac of vegetation not covered by snow [-]
    integer(ik4), pointer, contiguous, dimension(:) :: frac_veg_nosno_alb
    integer(ik4) :: fp, p, c   ! indices
    real(rk8) :: ol   ! thickness of canopy layer covered by snow (m)
    real(rk8) :: fb   ! fraction of canopy layer covered by snow

    if ( doalb ) then

      ! Assign local pointers to derived type scalar members (column-level)

      frac_sno => clm3%g%l%c%cps%frac_sno
      snow_depth  => clm3%g%l%c%cps%snow_depth

      ! Assign local pointers to derived type scalar members (pftlevel)

      pcolumn => clm3%g%l%c%p%column
      tlai    => clm3%g%l%c%p%pps%tlai
      tsai    => clm3%g%l%c%p%pps%tsai
      elai    => clm3%g%l%c%p%pps%elai
      esai    => clm3%g%l%c%p%pps%esai
      htop    => clm3%g%l%c%p%pps%htop
      hbot    => clm3%g%l%c%p%pps%hbot
      frac_veg_nosno_alb => clm3%g%l%c%p%pps%frac_veg_nosno_alb
      ivt     => clm3%g%l%c%p%itype

      do fp = 1, num_nolakep
        p = filter_nolakep(fp)
        c = pcolumn(p)

        ! need to update elai and esai only every albedo time step so do not
        ! have any inconsistency in lai and sai between SurfaceAlbedo calls
        ! (i.e., if albedos are not done every time step).
        ! leaf phenology
        ! Set leaf and stem areas based on day of year
        ! Interpolate leaf area index, stem area index, and vegetation heights
        ! between two monthly
        ! The weights below (timwt(1) and timwt(2)) were obtained by a call to
        ! routine InterpMonthlyVeg in subroutine NCARlsm.
        !                 Field   Monthly Values
        !                -------------------------
        ! leaf area index LAI  <- mlai1 and mlai2
        ! leaf area index SAI  <- msai1 and msai2
        ! top height      HTOP <- mhvt1 and mhvt2
        ! bottom height   HBOT <- mhvb1 and mhvb2

        tlai(p) = timwt(1)*mlai2t(p,1) + timwt(2)*mlai2t(p,2)
        tsai(p) = timwt(1)*msai2t(p,1) + timwt(2)*msai2t(p,2)
        htop(p) = timwt(1)*mhvt2t(p,1) + timwt(2)*mhvt2t(p,2)
        hbot(p) = timwt(1)*mhvb2t(p,1) + timwt(2)*mhvb2t(p,2)

        ! adjust lai and sai for burying by snow. if exposed lai and sai
        ! are less than 0.05, set equal to zero to prevent numerical
        ! problems associated with very small lai and sai.

        ! snow burial fraction for short vegetation (e.g. grasses) as in
        ! Wang and Zeng, 2007.

        if ( ivt(p) > noveg .and. ivt(p) <= nbrdlf_dcd_brl_shrub ) then
          ol = min( max(snow_depth(c)-hbot(p), 0._rk8), htop(p)-hbot(p))
          fb = 1._rk8 - ol / max(1.e-6_rk8, htop(p)-hbot(p))
        else
          !depth of snow required for complete burial of grasses
          fb = 1._rk8 - max(min(snow_depth(c),0.2_rk8),0._rk8)/0.2_rk8 ! 0.2m is assumed
        end if

        fb = min(max(fb,0.0_rk8),1.0_rk8)

        ! area weight by snow covered fraction
        elai(p) = max(tlai(p)*(1.0_rk8 - frac_sno(c)) &
               +tlai(p)*fb*frac_sno(c), 0.0_rk8)
        esai(p) = max(tsai(p)*(1.0_rk8 - frac_sno(c)) &
               +tsai(p)*fb*frac_sno(c), 0.0_rk8)
        if ( elai(p) < 0.05_rk8 ) elai(p) = 0.0_rk8
        if ( esai(p) < 0.05_rk8 ) esai(p) = 0.0_rk8
        if ( elai(p) > 10.0_rk8 ) elai(p) = 10._rk8
        if ( esai(p) > 10.0_rk8 ) esai(p) = 10._rk8

        ! Fraction of vegetation free of snow

        if ( (elai(p) + esai(p)) >= 0.05_rk8 ) then
          frac_veg_nosno_alb(p) = 1
        else
          frac_veg_nosno_alb(p) = 0
        end if
      end do ! end of pft loop
    end if  !end of if-doalb block
  end subroutine EcosystemDyn

#endif
  !
  ! Determine if 2 new months of data are to be read.
  !
  subroutine interpMonthlyVeg()
    implicit none
    integer(ik4) :: kyr   ! year (0, ...) for nstep+1
    integer(ik4) :: kmo   ! month (1, ..., 12)
    integer(ik4) :: kda   ! day of month (1, ..., 31)
    integer(ik4) :: ksec  ! seconds into current date for nstep+1
    integer(ik4) :: yr1, yr2
    real(rk8) :: t, rs, f1, f2
    integer(ik4), dimension(2) :: months ! months to be interpolated (1 to 12)

    call curr_date(nextdate,kyr,kmo,kda,ksec)
    t = real((kda-1)*86400+ksec,rk8) / &
                real(86400*ndaypm(kyr,kmo,nextdate%calendar),rk8)
    if ( t > 0.5_rk8 ) then
      months(1) = kmo
      months(2) = kmo + 1
    else
      months(1) = kmo - 1
      months(2) = kmo
    end if
    yr1 = kyr
    yr2 = kyr
    if ( months(1) < 1 ) then
      months(1) = 12
      yr1 = yr1 - 1
    end if
    if ( months(2) > 12 ) then
      months(2) = 1
      yr2 = yr2 + 1
    end if
    if ( InterpMonths1 /= months(1) ) then
      if (myid == italk) then
        write(stdout,*) 'Attempting to read monthly vegetation data .....'
        write(stdout,*) 'At ',trim(tochar(nextdate))
        write(stdout,*) 'Month = ',kmo,' Day = ',kda
        write(stdout,*) 'Reading months : ',months(1),months(2)
      end if
      call readMonthlyVegetation (fsurdat, months)
      InterpMonths1 = months(1)
    end if
    f1 = (ndaypm(yr1,months(1),nextdate%calendar)*0.5_rk8)*86400.0_rk8
    f2 = (ndaypm(yr2,months(2),nextdate%calendar)*0.5_rk8)*86400.0_rk8
    rs = (f1+f2)
    if ( t > 0.5_rk8 ) then
      timwt(2) = (real((kda-1)*86400+ksec,rk8) - f1)/rs
    else
      timwt(2) = (real((kda-1)*86400+ksec,rk8) + f1)/rs
    end if
    timwt(1) = 1.0_rk8-timwt(2)
    if ( any(timwt < 0.0_rk8) .or. any(timwt > 1.0_rk8) ) then
      call fatal(__FILE__,__LINE__, &
        'Incorrect interpolation weights!')
    end if
  end subroutine interpMonthlyVeg
  !
  !-----------------------------------------------------------------------
  ! read 12 months of veg data for dry deposition
  !-----------------------------------------------------------------------
  !
  subroutine readAnnualVegetation ( )
    implicit none
    type(clm_filetype) :: ncid  ! netcdf id
    ! 12 months of monthly lai from input data set
    real(rk8), pointer, contiguous, dimension(:,:) :: annlai
    ! lai read from input files
    real(rk8), pointer, contiguous, dimension(:,:) :: mlai
    integer(ik4) :: ier ! error code
    integer(ik4) :: g, k, l, p, ivt ! indices
    integer(ik4) :: ni, nj, ns   ! indices
    integer(ik4) :: begp, endp    ! beg and end local p index
    integer(ik4) :: begg, endg    ! beg and end local g index
    character(len=32) :: subname = 'readAnnualVegetation'

    annlai    => clm3%g%l%c%p%pps%annlai

    ! Determine necessary indices

    call get_proc_bounds(begg=begg,endg=endg,begp=begp,endp=endp)

    allocate(mlai(begg:endg,0:numpft), stat=ier)
    if (ier /= 0) then
      write(stderr,*)subname, 'allocation error '
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    if (myid == italk) then
      write (stdout,*) 'Attempting to read annual vegetation data .....'
    end if

    call clm_openfile(trim(fsurdat),ncid)
    call clm_check_dims(ncid,ni,nj)
    ns = ni*nj

!FAB obsolete check ??
!    if ( ldomain%ns /= ns .or. ldomain%ni /= ni .or. ldomain%nj /= nj ) then
!      write(stderr,*)trim(subname), 'ldomain and input file do not match dims '
!      write(stderr,*)trim(subname), 'ldomain%ni,ni,= ',ldomain%ni,ni
!      write(stderr,*)trim(subname), 'ldomain%nj,nj,= ',ldomain%nj,nj
!      write(stderr,*)trim(subname), 'ldomain%ns,ns,= ',ldomain%ns,ns
!      call fatal(__FILE__,__LINE__,'clm now stopping')
!    end if
    if ( .not. clm_check_dimlen(ncid,'lsmpft',numpft+1) ) then
      call fatal(__FILE__,__LINE__, &
              trim(subname)//' ERROR: lsmpft not equal to numpft+1')
    end if

    do k = 1, 12   !! loop over months and read vegetated data
      call clm_readvar(ncid,'MONTHLY_LAI',mlai,gcomm_gridcell,k)
      !! store data directly in clmtype structure
      !! only vegetated pfts have nonzero values
      !! Assign lai/sai/hgtt/hgtb to the top [maxpatch_pft] pfts
      !! as determined in subroutine surfrd
      do p = begp, endp
        g = clm3%g%l%c%p%gridcell(p)
        ivt = clm3%g%l%c%p%itype(p)
        if ( ivt /= noveg ) then     !! vegetated pft
          do l = 0, numpft
            if ( l == ivt ) then
              annlai(k,p) = mlai(g,l)
            end if
          end do
        else                       !! non-vegetated pft
          annlai(k,p) = 0._rk8
        end if
      end do   ! end of loop over pfts
    end do ! months loop
    call clm_closefile(ncid)
    deallocate(mlai)
  end subroutine readAnnualVegetation
  !
  ! Read monthly vegetation data for two consec. months.
  !
  subroutine readMonthlyVegetation (fveg, months)
    implicit none
    ! file with monthly vegetation data
    character(len=*), intent(in) :: fveg
    ! months to be interpolated (1 to 12)
    integer(ik4), dimension(2), intent(in) :: months
    type(clm_filetype) :: ncid            ! netcdf id
    integer(ik4) :: g, k, l, p, ivt ! indices
    integer(ik4) :: begp, endp    ! beg and end local p index
    integer(ik4) :: begg, endg    ! beg and end local g index
    integer(ik4) :: ier                        ! error code
    ! lai read from input files
    real(rk8), pointer, contiguous, dimension(:,:) :: mlai
    ! sai read from input files
    real(rk8), pointer, contiguous, dimension(:,:) :: msai
    ! top vegetation height
    real(rk8), pointer, contiguous, dimension(:,:) :: mhgtt
    ! bottom vegetation height
    real(rk8), pointer, contiguous, dimension(:,:) :: mhgtb
    ! difference between lai month one and month two
    real(rk8), pointer, contiguous, dimension(:) :: mlaidiff
    character(len=32) :: subname = 'readMonthlyVegetation'

    ! Determine necessary indices

    call get_proc_bounds(begg=begg,endg=endg,begp=begp,endp=endp)

    allocate(mlai(begg:endg,0:numpft), &
             msai(begg:endg,0:numpft), &
             mhgtt(begg:endg,0:numpft), &
             mhgtb(begg:endg,0:numpft), &
             stat=ier)
    if (ier /= 0) then
      write(stderr,*)subname, 'allocation big error '
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! ----------------------------------------------------------------------
    ! Open monthly vegetation file
    ! Read data and convert from gridcell to pft data
    ! ----------------------------------------------------------------------

    call clm_openfile(fveg,ncid)

    do k = 1, 2   !loop over months and read vegetated data
      call clm_readvar(ncid,'MONTHLY_LAI',mlai,gcomm_gridcell,months(k))
      call clm_readvar(ncid,'MONTHLY_SAI',msai,gcomm_gridcell,months(k))
      call clm_readvar(ncid,'MONTHLY_HEIGHT_TOP',mhgtt,gcomm_gridcell,months(k))
      call clm_readvar(ncid,'MONTHLY_HEIGHT_BOT',mhgtb,gcomm_gridcell,months(k))
      ! Store data directly in clmtype structure
      ! only vegetated pfts have nonzero values
      ! Assign lai/sai/hgtt/hgtb to the top [maxpatch_pft] pfts
      ! as determined in subroutine surfrd
      do p = begp, endp
        g = clm3%g%l%c%p%gridcell(p)
        ivt = clm3%g%l%c%p%itype(p)
        if ( ivt /= noveg ) then     ! vegetated pft
          do l = 0, numpft
            if ( l == ivt ) then
              mlai2t(p,k) = mlai(g,l)
              msai2t(p,k) = msai(g,l)
              mhvt2t(p,k) = mhgtt(g,l)
              mhvb2t(p,k) = mhgtb(g,l)
            end if
          end do
        else                        ! non-vegetated pft
          mlai2t(p,k) = 0._rk8
          msai2t(p,k) = 0._rk8
          mhvt2t(p,k) = 0._rk8
          mhvb2t(p,k) = 0._rk8
        end if
      end do   ! end of loop over pfts
    end do   ! end of loop over months
    call clm_closefile(ncid)
    if (myid == italk) then
      k = 2
      write(stdout,*) 'Successfully read monthly vegetation data for'
      write(stdout,*) 'month ', months(k)
    end if
    deallocate(mlai, msai, mhgtt, mhgtb)
    mlaidiff => clm3%g%l%c%p%pps%mlaidiff
    do p = begp, endp
      mlaidiff(p) = mlai2t(p,1)-mlai2t(p,2)
    enddo
  end subroutine readMonthlyVegetation

end module mod_clm_staticecosysdyn
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
