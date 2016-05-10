module mod_clm_dust
  !
  ! Routines in this module calculate Dust mobilization and dry deposition
  ! for dust.
  ! Simulates dust mobilization due to wind from the surface into the
  ! lowest atmospheric layer. On output flx_mss_vrt_dst(ndst) is the surface
  ! dust emission (kg/m**2/s) [ + = to atm].
  ! Calculates the turbulent component of dust dry deposition, (the turbulent
  ! deposition velocity through the lowest atmospheric layer).
  ! CAM will calculate the settling velocity through the whole atmospheric
  ! column. The two calculations will determine the dust dry deposition flux
  ! to the surface.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_stdio
  use mod_clm_type
  use mod_clm_varpar  , only : dst_src_nbr, ndst, sz_nbr
  use mod_clm_varcon  , only : grav, istsoil
  use mod_clm_varcon  , only : istcrop
  use mod_clm_subgridave, only: p2l
  use mod_clm_varcon , only: spval , denh2o
  use mod_clm_atmlnd   , only : clm_a2l

  implicit none

  private

  save
  public :: Dustini        ! Initialize variables used in subroutine Dust
  public :: DustEmission   ! Dust mobilization
  public :: DustDryDep     ! Turbulent dry deposition for dust

  real(rkx) :: ovr_src_snk_mss(dst_src_nbr,ndst)
  !Factor in saltation computation (named as in Charlie's code)
  real(rkx) :: tmp1
  real(rkx) :: dmt_vwr(ndst) ![m] Mass-weighted mean diameter resolved
  real(rkx) :: stk_crc(ndst) ![frc] Correction to Stokes settling velocity
  real(rkx) :: dns_aer       ![kg m-3] Aerosol density

  contains
  !
  ! Dust mobilization. This code simulates dust mobilization due to wind
  ! from the surface into the lowest atmospheric layer
  ! On output flx_mss_vrt_dst(ndst) is the surface dust emission
  ! (kg/m**2/s) [ + = to atm]
  ! Source: C. Zender's dust model
  !
  subroutine DustEmission(lbp,ubp,lbl,ubl,num_nolakep,filter_nolakep)
    implicit none
    integer(ik4), intent(in) :: lbp,ubp,ubl,lbl ! pft bounds
    ! number of column non-lake points in pft filter
    integer(ik4), intent(in) :: num_nolakep
    ! pft filter for non-lake points
    integer(ik4), intent(in) :: filter_nolakep(num_nolakep)
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer :: pactive(:)
    integer(ik4) , pointer :: pcolumn(:)   ! pft's column index
    integer(ik4) , pointer :: plandunit(:) ! pft's landunit index
    integer(ik4) , pointer :: pgridcell(:) ! pft's gridcell index
    integer(ik4) , pointer :: ityplun(:)   ! landunit type
    ! one-sided leaf area index, no burying by snow
    real(rkx), pointer :: tlai(:)
    ! one-sided stem area index, no burying by snow
    real(rkx), pointer :: tsai(:)
    ! fraction of ground covered by snow (0 to 1)
    real(rkx), pointer :: frac_sno(:)
    ! threshold gravimetric soil moisture based on clay content
    real(rkx), pointer :: gwc_thr(:)
    real(rkx), pointer :: forc_rho(:)   ! density (kg/m**3)
    real(rkx), pointer :: fv(:) ! friction velocity (m/s) (for dust model)
    real(rkx), pointer :: u10(:) ! 10-m wind (m/s) (created for dust model)
    real(rkx), pointer :: mbl_bsn_fct(:)     ! basin factor
    ! [frc] Mass fraction clay limited to 0.20
    real(rkx), pointer :: mss_frc_cly_vld(:)
    ! volumetric soil water (0<=h2osoi_vol<=watsat)
    real(rkx), pointer :: h2osoi_vol(:,:)
    real(rkx), pointer :: h2osoi_liq(:,:)  ! liquid soil water (kg/m2)
    real(rkx), pointer :: h2osoi_ice(:,:)  ! frozen soil water (kg/m2)
    real(rkx), pointer :: watsat(:,:)      ! saturated volumetric soil water

    ! surface dust emission (kg/m**2/s)
    real(rkx), pointer :: flx_mss_vrt_dst(:,:)
    ! total dust flux into atmosphere
    real(rkx), pointer :: flx_mss_vrt_dst_tot(:)

    integer(ik4)  :: fp,p,c,l,g,m,n  ! indices
    ! fraction of total water that is liquid
    real(rkx) :: liqfrac
    ! [frc] Wind friction threshold over wind friction
    real(rkx) :: wnd_frc_rat
    ! [m s-1] Friction velocity increase from saltatn
    real(rkx) :: wnd_frc_slt_dlt
    ! [m s-1] Reference windspeed excess over threshld
    real(rkx) :: wnd_rfr_dlt
    real(rkx) :: dst_slt_flx_rat_ttl
    real(rkx) :: flx_mss_hrz_slt_ttl
    real(rkx) :: flx_mss_vrt_dst_ttl(lbp:ubp)
    real(rkx) :: frc_thr_wet_fct
    real(rkx) :: frc_thr_rgh_fct
    real(rkx) :: wnd_frc_thr_slt
    real(rkx) :: wnd_rfr_thr_slt
    real(rkx) :: wnd_frc_slt
    real(rkx) :: lnd_frc_mbl(lbp:ubp)
    real(rkx) :: bd
    real(rkx) :: gwc_sfc
    real(rkx) :: ttlai(lbp:ubp)
    real(rkx) :: tlai_lu(lbl:ubl)

    ! [frc] Saltation constant
    real(rkx), parameter :: cst_slt = 2.61_rkx
    ! [frc] Empir. mass flx tuning eflx_lh_vegt
    real(rkx), parameter :: flx_mss_fdg_fct = 5.0e-4_rkx
    ! [m2 m-2] VAI threshold quenching dust mobilization
    real(rkx), parameter :: vai_mbl_thr = 0.3_rkx
    real(rkx), pointer :: wtlunit(:)  ! weight of pft relative to landunit
    real(rkx) :: sumwt(lbl:ubl)       ! sum of weights
    logical  :: found                 ! temporary for error check
    integer(ik4)  :: index

    ! Assign local pointers to derived type scalar members (gridcell-level)

    forc_rho        => clm_a2l%forc_rho

    ! Assign local pointers to derived type scalar members (landunit-level)

    ityplun         => clm3%g%l%itype

    ! Assign local pointers to derived type scalar members (column-level)

    frac_sno        => clm3%g%l%c%cps%frac_sno
    gwc_thr         => clm3%g%l%c%cps%gwc_thr
    mbl_bsn_fct     => clm3%g%l%c%cps%mbl_bsn_fct
    mss_frc_cly_vld => clm3%g%l%c%cps%mss_frc_cly_vld
    h2osoi_vol      => clm3%g%l%c%cws%h2osoi_vol
    h2osoi_liq      => clm3%g%l%c%cws%h2osoi_liq
    h2osoi_ice      => clm3%g%l%c%cws%h2osoi_ice
    watsat          => clm3%g%l%c%cps%watsat

    ! Assign local pointers to derived type scalar members (pft-level)

    pactive         => clm3%g%l%c%p%active
    pgridcell       => clm3%g%l%c%p%gridcell
    plandunit       => clm3%g%l%c%p%landunit
    pcolumn         => clm3%g%l%c%p%column
    tlai            => clm3%g%l%c%p%pps%tlai
    tsai            => clm3%g%l%c%p%pps%tsai
    fv              => clm3%g%l%c%p%pps%fv
    u10             => clm3%g%l%c%p%pps%u10
    flx_mss_vrt_dst => clm3%g%l%c%p%pdf%flx_mss_vrt_dst
    flx_mss_vrt_dst_tot => clm3%g%l%c%p%pdf%flx_mss_vrt_dst_tot
   !local pointers from subgridAveMod/p2l_1d
    wtlunit         => clm3%g%l%c%p%wtlunit

    ttlai(:) = 0._rkx
    ! make lai average at landunit level
    do fp = 1,num_nolakep
      p = filter_nolakep(fp)
      ttlai(p) = tlai(p)+tsai(p)
    enddo

    tlai_lu(:) = spval
    sumwt(:) = 0._rkx
    do p = lbp,ubp
      if (ttlai(p) /= spval .and. pactive(p) .and. wtlunit(p) /= 0._rkx) then
        c = pcolumn(p)
        l = plandunit(p)
        if (sumwt(l) == 0._rkx) tlai_lu(l) = 0._rkx
        tlai_lu(l) = tlai_lu(l) + ttlai(p) * wtlunit(p)
        sumwt(l) = sumwt(l) + wtlunit(p)
      end if
    end do
    found = .false.
    do l = lbl,ubl
      if (sumwt(l) > 1.0_rkx + 1.e-4_rkx) then
        found = .true.
        index = l
        exit
      else if (sumwt(l) /= 0._rkx) then
        tlai_lu(l) = tlai_lu(l)/sumwt(l)
      end if
    end do
    if (found) then
      write(stderr,*) 'p2l_1d error: sumwt is greater than 1.0 at l= ',index
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! Loop through pfts

    ! initialize variables which get passed to the atmosphere
    flx_mss_vrt_dst(lbp:ubp,:)=0._rkx

    do fp = 1,num_nolakep
      p = filter_nolakep(fp)
      c = pcolumn(p)
      l = plandunit(p)

      ! the following code from subr. lnd_frc_mbl_get was adapted for lsm use
      ! purpose: return fraction of each gridcell suitable for dust mobilization
      ! the "bare ground" fraction of the current sub-gridscale cell decreases
      ! linearly from 1 to 0 as VAI(=tlai+tsai) increases from 0 to vai_mbl_thr
      ! if ice sheet, wetland, or lake, no dust allowed

      if (ityplun(l) == istsoil .or. ityplun(l) == istcrop) then
        if (tlai_lu(l) < vai_mbl_thr) then
          lnd_frc_mbl(p) = 1.0_rkx - (tlai_lu(l))/vai_mbl_thr
        else
          lnd_frc_mbl(p) = 0.0_rkx
        endif
        lnd_frc_mbl(p) = lnd_frc_mbl(p) * (1.0_rkx - frac_sno(c))
      else
        lnd_frc_mbl(p) = 0.0_rkx
      end if
    end do

    do fp = 1,num_nolakep
      p = filter_nolakep(fp)
      if (lnd_frc_mbl(p)>1.0_rkx .or. lnd_frc_mbl(p)<0.0_rkx) then
        write(stderr,*) &
                'Error dstmbl: pft= ',p,' lnd_frc_mbl(p)= ',lnd_frc_mbl(p)
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
    end do

    ! reset history output variables before next if-statement
    ! to avoid output = inf

    do fp = 1,num_nolakep
      p = filter_nolakep(fp)
      flx_mss_vrt_dst_tot(p) = 0.0_rkx
    end do
    do n = 1, ndst
      do fp = 1,num_nolakep
        p = filter_nolakep(fp)
        flx_mss_vrt_dst(p,n) = 0.0_rkx
      end do
    end do

    do fp = 1,num_nolakep
      p = filter_nolakep(fp)
      c = pcolumn(p)
      l = plandunit(p)
      g = pgridcell(p)

      ! only perform the following calculations if lnd_frc_mbl is non-zero

      if (lnd_frc_mbl(p) > 0.0_rkx) then

        ! the following comes from subr. frc_thr_rgh_fct_get
        ! purpose: compute factor by which surface roughness increases threshold
        !          friction velocity (currently a constant)

        frc_thr_rgh_fct = 1.0_rkx

        ! the following comes from subr. frc_thr_wet_fct_get
        ! purpose: compute factor by which soil moisture increases
        ! threshold friction velocity
        ! adjust threshold velocity for inhibition by moisture
        ! modified 4/5/2002 (slevis) to use gravimetric instead of volumetric
        ! water content

        ![kg m-3] Bulk density of dry surface soil
        bd = (1._rkx-watsat(c,1))*2.7e3_rkx
        ![kg kg-1] Gravimetric H2O cont
        gwc_sfc = h2osoi_vol(c,1)*denh2o/bd
        if (gwc_sfc > gwc_thr(c)) then
          frc_thr_wet_fct = sqrt(1.0_rkx + &
                  1.21_rkx * (100.0_rkx*(gwc_sfc - gwc_thr(c)))**0.68_rkx)
        else
          frc_thr_wet_fct = 1.0_rkx
        end if

        ! slevis: adding liqfrac here, because related to effects
        ! from soil water

        liqfrac = max( 0.0_rkx, min( 1.0_rkx,  &
                h2osoi_liq(c,1) / (h2osoi_ice(c,1)+h2osoi_liq(c,1)+1.0e-6_rkx) ) )

        ! the following lines come from subr. dst_mbl
        ! purpose: adjust threshold friction velocity to acct for moisture and
        !          roughness. The ratio tmp1 / sqrt(forc_rho) comes from
        !          subr. wnd_frc_thr_slt_get which computes dry threshold
        !          friction velocity for saltation

        wnd_frc_thr_slt = tmp1 / sqrt(forc_rho(g)) * &
                frc_thr_wet_fct * frc_thr_rgh_fct

        ! reset these variables which will be updated in the following if-block

        wnd_frc_slt = fv(p)
        flx_mss_hrz_slt_ttl = 0.0_rkx
        flx_mss_vrt_dst_ttl(p) = 0.0_rkx

        ! the following line comes from subr. dst_mbl
        ! purpose: threshold saltation wind speed

        wnd_rfr_thr_slt = u10(p) * wnd_frc_thr_slt / fv(p)

        ! the following if-block comes from subr. wnd_frc_slt_get
        ! purpose: compute the saltating friction velocity
        ! theory: saltation roughens the boundary layer, AKA "Owen's effect"

        if (u10(p) >= wnd_rfr_thr_slt) then
          wnd_rfr_dlt = u10(p) - wnd_rfr_thr_slt
          wnd_frc_slt_dlt = 0.003_rkx * wnd_rfr_dlt * wnd_rfr_dlt
          wnd_frc_slt = fv(p) + wnd_frc_slt_dlt
        end if

        ! the following comes from subr. flx_mss_hrz_slt_ttl_Whi79_get
        ! purpose: compute vertically integrated streamwise
        ! mass flux of particles

        if (wnd_frc_slt > wnd_frc_thr_slt) then
          wnd_frc_rat = wnd_frc_thr_slt / wnd_frc_slt
          flx_mss_hrz_slt_ttl = cst_slt * forc_rho(g) * (wnd_frc_slt**3.0_rkx) * &
               (1.0_rkx - wnd_frc_rat) * (1.0_rkx + wnd_frc_rat) * &
               (1.0_rkx + wnd_frc_rat) / grav

          ! the following loop originates from subr. dst_mbl
          ! purpose: apply land sfc and veg limitations and global tuning factor
          ! slevis: multiply flx_mss_hrz_slt_ttl by liqfrac to incude the effect
          ! of frozen soil

          flx_mss_hrz_slt_ttl = flx_mss_hrz_slt_ttl * &
                  lnd_frc_mbl(p) * mbl_bsn_fct(c) * &
                  flx_mss_fdg_fct * liqfrac
        end if

        ! the following comes from subr. flx_mss_vrt_dst_ttl_MaB95_get
        ! purpose: diagnose total vertical mass flux of dust from vertically
        !          integrated streamwise mass flux

        dst_slt_flx_rat_ttl = 100.0_rkx * &
                exp( log(10.0_rkx) * (13.4_rkx * mss_frc_cly_vld(c) - 6.0_rkx) )
        flx_mss_vrt_dst_ttl(p) = flx_mss_hrz_slt_ttl * dst_slt_flx_rat_ttl

      end if   ! lnd_frc_mbl > 0.0
    end do

    ! the following comes from subr. flx_mss_vrt_dst_prt in C. Zender's code
    ! purpose: partition total vertical mass flux of dust into transport bins

    do n = 1, ndst
      do m = 1, dst_src_nbr
        do fp = 1,num_nolakep
          p = filter_nolakep(fp)
          if (lnd_frc_mbl(p) > 0.0_rkx) then
            flx_mss_vrt_dst(p,n) = flx_mss_vrt_dst(p,n) + &
                    ovr_src_snk_mss(m,n) * flx_mss_vrt_dst_ttl(p)
          end if
        end do
      end do
    end do

    do n = 1, ndst
      do fp = 1,num_nolakep
        p = filter_nolakep(fp)
        if (lnd_frc_mbl(p) > 0.0_rkx) then
          flx_mss_vrt_dst_tot(p) = flx_mss_vrt_dst_tot(p) + flx_mss_vrt_dst(p,n)
        end if
      end do
    end do
  end subroutine DustEmission
  !
  ! Determine Turbulent dry deposition for dust. Calculate the turbulent
  ! component of dust dry deposition, (the turbulent deposition velocity
  ! through the lowest atmospheric layer. CAM will calculate the settling
  ! velocity through the whole atmospheric column. The two calculations
  ! will determine the dust dry deposition flux to the surface.
  ! Note: Same process should occur over oceans. For the coupled CESM,
  ! we may find it more efficient to let CAM calculate the turbulent dep
  ! velocity over all surfaces. This would require passing the
  ! aerodynamic resistance, ram(1), and the friction velocity, fv, from
  ! the land to the atmosphere component. In that case, dustini need not
  ! calculate particle diamter (dmt_vwr) and particle density (dns_aer).
  ! Source: C. Zender's dry deposition code
  !
  subroutine DustDryDep (lbp, ubp)
    use mod_constants , only : boltzk
    use mod_clm_varcon , only : rpi , rair
    use mod_clm_atmlnd , only : clm_a2l
    implicit none
    integer(ik4), intent(in) :: lbp, ubp     ! pft bounds
    ! true=>do computations on this pft (see reweightMod for details)
    logical , pointer :: pactive(:)
    integer(ik4) , pointer :: pgridcell(:)   ! pft's gridcell index
    real(rkx), pointer :: forc_t(:)      ! atm temperature (K)
    real(rkx), pointer :: forc_pbot(:)   ! atm pressure (Pa)
    real(rkx), pointer :: forc_rho(:)    ! atm density (kg/m**3)
    real(rkx), pointer :: fv(:)          ! friction velocity (m/s)
    real(rkx), pointer :: ram1(:)        ! aerodynamical resistance (s/m)
    real(rkx), pointer :: vlc_trb(:,:)   ! Turbulent deposn velocity (m/s)
    real(rkx), pointer :: vlc_trb_1(:)   ! Turbulent deposition velocity 1
    real(rkx), pointer :: vlc_trb_2(:)   ! Turbulent deposition velocity 2
    real(rkx), pointer :: vlc_trb_3(:)   ! Turbulent deposition velocity 3
    real(rkx), pointer :: vlc_trb_4(:)   ! Turbulent deposition velocity 4
    integer(ik4)  :: p,g,m                 ! indices
    real(rkx) :: vsc_dyn_atm(lbp:ubp)  ! [kg m-1 s-1] Dynamic viscosity of air
    ! [m2 s-1] Kinematic viscosity of atmosphere
    real(rkx) :: vsc_knm_atm(lbp:ubp)
    ! [frc] Sfc-dep exponent for aerosol-diffusion dependence on Schmidt number
    real(rkx) :: shm_nbr_xpn
    real(rkx) :: shm_nbr      ! [frc] Schmidt number
    real(rkx) :: stk_nbr      ! [frc] Stokes number
    real(rkx) :: mfp_atm      ! [m] Mean free path of air
    real(rkx) :: dff_aer      ! [m2 s-1] Brownian diffusivity of particle
    real(rkx) :: rss_trb      ! [s m-1] Resistance to turbulent deposition
    real(rkx) :: slp_crc(lbp:ubp,ndst) ! [frc] Slip correction factor
    real(rkx) :: vlc_grv(lbp:ubp,ndst) ! [m s-1] Settling velocity
    real(rkx) :: rss_lmn(lbp:ubp,ndst) ! [s m-1] Quasi-laminar layer resistance
    real(rkx) :: tmp                   ! temporary

    real(rkx),parameter::shm_nbr_xpn_lnd=-2._rkx/3._rkx ![frc] shm_nbr_xpn over land

    ! Assign local pointers to derived type members (gridcell-level)

    forc_pbot => clm_a2l%forc_pbot
    forc_rho  => clm_a2l%forc_rho
    forc_t    => clm_a2l%forc_t

    ! Assign local pointers to derived type members (pft-level)

    pactive   => clm3%g%l%c%p%active
    pgridcell => clm3%g%l%c%p%gridcell
    fv        => clm3%g%l%c%p%pps%fv
    ram1      => clm3%g%l%c%p%pps%ram1
    vlc_trb   => clm3%g%l%c%p%pdf%vlc_trb
    vlc_trb_1 => clm3%g%l%c%p%pdf%vlc_trb_1
    vlc_trb_2 => clm3%g%l%c%p%pdf%vlc_trb_2
    vlc_trb_3 => clm3%g%l%c%p%pdf%vlc_trb_3
    vlc_trb_4 => clm3%g%l%c%p%pdf%vlc_trb_4

    do p = lbp,ubp
      if (pactive(p)) then
        g = pgridcell(p)

        ! from subroutine dst_dps_dry
        ! (consider adding sanity checks from line 212)
        ! when code asks to use midlayer density, pressure, temperature,
        ! I use the data coming in from the atmosphere, ie
        ! forc_t, forc_pbot, forc_rho

        ! Quasi-laminar layer resistance: call rss_lmn_get
        ! Size-independent thermokinetic properties

        vsc_dyn_atm(p) = 1.72e-5_rkx * ((forc_t(g)/273.0_rkx)**1.5_rkx) * 393.0_rkx / &
               (forc_t(g)+120.0_rkx)      ![kg m-1 s-1] RoY94 p. 102
        mfp_atm = 2.0_rkx * vsc_dyn_atm(p) / &   ![m] SeP97 p. 455
               (forc_pbot(g)*sqrt(8.0_rkx/(rpi*rair*forc_t(g))))
        ![m2 s-1] Kinematic viscosity of air
        vsc_knm_atm(p) = vsc_dyn_atm(p) / forc_rho(g)

        do m = 1, ndst
          ![frc] Slip correction factor SeP97 p. 464
          if ( (1.1_rkx*dmt_vwr(m)/(2.0_rkx*mfp_atm)) > 25.0 ) then
            slp_crc(p,m) = 1.0_rkx + 2.0_rkx * mfp_atm * &
                    1.257_rkx / dmt_vwr(m)
          else
            slp_crc(p,m) = 1.0_rkx + 2.0_rkx * mfp_atm * &
                    (1.257_rkx+0.4_rkx*exp(-1.1_rkx*dmt_vwr(m)/(2.0_rkx*mfp_atm))) / &
                    dmt_vwr(m)
          end if
          ![m s-1] Stokes' settling velocity SeP97 p. 466
          vlc_grv(p,m) = (1.0_rkx/18.0_rkx) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
                  grav * slp_crc(p,m) / vsc_dyn_atm(p)
          ![m s-1] Correction to Stokes settling velocity
          vlc_grv(p,m) = vlc_grv(p,m) * stk_crc(m)
        end do
      end if
    end do

    do m = 1, ndst
      do p = lbp,ubp
        if (pactive(p)) then
          g = pgridcell(p)

          ![frc] SeP97 p.965
          stk_nbr = vlc_grv(p,m) * fv(p) * fv(p) / (grav * vsc_knm_atm(p))
          dff_aer = boltzk * forc_t(g) * slp_crc(p,m) / &    ![m2 s-1]
              (3.0_rkx*rpi * vsc_dyn_atm(p) * dmt_vwr(m))      !SeP97 p.474
          shm_nbr = vsc_knm_atm(p) / dff_aer                 ![frc] SeP97 p.972
             shm_nbr_xpn = shm_nbr_xpn_lnd                   ![frc]

          ! fxm: Turning this on dramatically reduces
          ! deposition velocity in low wind regimes
          ! Schmidt number exponent is -2/3 over solid surfaces and
          ! -1/2 over liquid surfaces SlS80 p. 1014
          ! if (oro(i)==0.0) shm_nbr_xpn=shm_nbr_xpn_ocn else
          ! shm_nbr_xpn=shm_nbr_xpn_lnd
          ! [frc] Surface-dependent exponent for aerosol-diffusion
          ! dependence on Schmidt #

          if ( 3.0_rkx/stk_nbr < 25 ) then
            tmp = shm_nbr**shm_nbr_xpn + 10.0_rkx**(-3.0_rkx/stk_nbr)
          else
            tmp = shm_nbr**shm_nbr_xpn
          end if
          rss_lmn(p,m) = 1.0_rkx / (tmp * fv(p)) ![s m-1] SeP97 p.972,965
        end if
      end do
    end do

    ! Lowest layer: Turbulent deposition (CAM will calc. gravitational dep)

    do m = 1, ndst
      do p = lbp,ubp
        if (pactive(p)) then
          rss_trb = ram1(p) + rss_lmn(p,m) + &
                  ram1(p) * rss_lmn(p,m) * vlc_grv(p,m) ![s m-1]
          vlc_trb(p,m) = 1.0_rkx / rss_trb  ![m s-1]
        end if
      end do
    end do

    do p = lbp,ubp
      if (pactive(p)) then
        vlc_trb_1(p) = vlc_trb(p,1)
        vlc_trb_2(p) = vlc_trb(p,2)
        vlc_trb_3(p) = vlc_trb(p,3)
        vlc_trb_4(p) = vlc_trb(p,4)
      end if
    end do

  end subroutine DustDryDep
  !
  ! Compute source efficiency factor from topography
  ! Initialize other variables used in subroutine Dust:
  ! ovr_src_snk_mss(m,n) and tmp1.
  ! Define particle diameter and density needed by atm model
  ! as well as by dry dep model
  ! Source: Paul Ginoux (for source efficiency factor)
  ! Modifications by C. Zender and later by S. Levis
  ! Rest of subroutine from C. Zender's dust model
  !
  subroutine Dustini()
    use mod_clm_varcon , only : rpi , rair
    use mod_clm_decomp , only : get_proc_bounds
    implicit none
    real(rkx), pointer :: mbl_bsn_fct(:) !basin factor
    integer(ik4)  :: c,l,m,n              ! indices
    real(rkx) :: ovr_src_snk_frc
    real(rkx) :: sqrt2lngsdi             ! [frc] Factor in erf argument
    real(rkx) :: lndmaxjovrdmdni         ! [frc] Factor in erf argument
    real(rkx) :: lndminjovrdmdni         ! [frc] Factor in erf argument
    ! [frc] Threshold friction Reynolds number approximation for optimal size
    real(rkx) :: ryn_nbr_frc_thr_prx_opt
    ! [frc] Threshold friction Reynolds factor for saltation calculation
    real(rkx) :: ryn_nbr_frc_thr_opt_fnc
    ! Interpartical cohesive forces factor for saltation calc
    real(rkx) :: icf_fct
    ! Density ratio factor for saltation calculation
    real(rkx) :: dns_fct
    real(rkx) :: dmt_min(ndst)  ! [m] Size grid minimum
    real(rkx) :: dmt_max(ndst)  ! [m] Size grid maximum
    real(rkx) :: dmt_ctr(ndst)  ! [m] Diameter at bin center
    real(rkx) :: dmt_dlt(ndst)  ! [m] Width of size bin
    real(rkx) :: slp_crc(ndst)  ! [frc] Slip correction factor
    real(rkx) :: vlm_rsl(ndst)  ! [m3 m-3] Volume concentration resolved
    real(rkx) :: vlc_stk(ndst)  ! [m s-1] Stokes settling velocity
    real(rkx) :: vlc_grv(ndst)  ! [m s-1] Settling velocity
    real(rkx) :: ryn_nbr_grv(ndst) ! [frc] Reynolds number at terminal velocity
    real(rkx) :: cff_drg_grv(ndst) ! [frc] Drag coefficient at terminal velocity
    real(rkx) :: tmp      ! temporary
    real(rkx) :: ln_gsd   ! [frc] ln(gsd)
    real(rkx) :: gsd_anl  ! [frc] Geometric standard deviation
    real(rkx) :: dmt_vma  ! [m] Mass median diameter analytic She84 p.75 Tabl.1
    real(rkx) :: dmt_nma  ! [m] Number median particle diameter
    real(rkx) :: lgn_dst  ! Lognormal distribution at sz_ctr
    real(rkx) :: eps_max  ! [frc] Relative accuracy for convergence
    real(rkx) :: eps_crr  ! [frc] Current relative accuracy
    real(rkx) :: itr_idx  ! [idx] Counting index
    real(rkx) :: dns_mdp  ! [kg m-3] Midlayer density
    real(rkx) :: mfp_atm  ! [m] Mean free path of air
    real(rkx) :: vsc_dyn_atm  ! [kg m-1 s-1] Dynamic viscosity of air
    real(rkx) :: vsc_knm_atm  ! [kg m-1 s-1] Kinematic viscosity of air
    real(rkx) :: vlc_grv_old  ! [m s-1] Previous gravitational settling velocity
    real(rkx) :: series_ratio ! Factor for logarithmic grid
    real(rkx) :: lngsdsqrttwopi_rcp   ! Factor in lognormal distribution
    real(rkx) :: sz_min(sz_nbr)       ! [m] Size Bin minima
    real(rkx) :: sz_max(sz_nbr)       ! [m] Size Bin maxima
    real(rkx) :: sz_ctr(sz_nbr)       ! [m] Size Bin centers
    real(rkx) :: sz_dlt(sz_nbr)       ! [m] Size Bin widths

    ! constants
    ! [m] Mass median diameter
    real(rkx) :: dmt_vma_src(dst_src_nbr) =    &
         (/ 0.832e-6_rkx , 4.82e-6_rkx , 19.38e-6_rkx /)        ! BSM96 p. 73 Table 2
    ! [frc] Geometric std deviation
    real(rkx) :: gsd_anl_src(dst_src_nbr) =    &
         (/ 2.10_rkx     ,  1.90_rkx   , 1.60_rkx     /)        ! BSM96 p. 73 Table 2
    ! [frc] Mass fraction
    real(rkx) :: mss_frc_src(dst_src_nbr) =    &
         (/ 0.036_rkx, 0.957_rkx, 0.007_rkx /)                  ! BSM96 p. 73 Table 2
    ! [m] Particle diameter grid
    real(rkx) :: dmt_grd(5) =                  &
         (/ 0.1e-6_rkx, 1.0e-6_rkx, 2.5e-6_rkx, 5.0e-6_rkx, 10.0e-6_rkx /)
    ! [m] Optim diam for saltation
    real(rkx), parameter :: dmt_slt_opt = 75.0e-6_rkx
    ! [kg m-3] Density of optimal saltation particles
    real(rkx), parameter :: dns_slt = 2650.0_rkx

    ! declare erf intrinsic function
    real(rkx) :: dum     ! dummy variable for erf test
    integer(ik4) :: begp, endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc, endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl, endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg, endg ! per-proc gridcell ending gridcell indices

#ifdef __PGI
    real(rkx) , external :: erf
#endif
    ! Assign local pointers to derived type scalar members (column-level)

    mbl_bsn_fct => clm3%g%l%c%cps%mbl_bsn_fct

    ! Sanity check on erf: erf() in SGI /usr/lib64/mips4/libftn.so is bogus

    dum = 1.0_rkx
    if (abs(0.8427_rkx-erf(dum))/0.8427_rkx>0.001_rkx) then
       write(stderr,*) 'erf(1.0) = ',erf(dum)
       write(stderr,*) 'Dustini: Error function error'
       call fatal(__FILE__,__LINE__,'clm now stopping')
    end if
    dum = 0.0_rkx
    if (erf(dum) /= 0.0_rkx) then
       write(stderr,*) 'erf(0.0) = ',erf(dum)
       write(stderr,*) 'Dustini: Error function error'
       call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! the following comes from (1) szdstlgn.F subroutine ovr_src_snk_frc_get
    !                      and (2) dstszdst.F subroutine dst_szdst_ini
    ! purpose(1): given one set (the "source") of lognormal distributions,
    !             and one set of bin boundaries (the "sink"), compute and return
    !             the overlap factors between the source and sink distributions
    ! purpose(2): set important statistics of size distributions

    do m = 1, dst_src_nbr
      sqrt2lngsdi = sqrt(2.0_rkx) * log(gsd_anl_src(m))
      do n = 1, ndst
        lndmaxjovrdmdni = log(dmt_grd(n+1)/dmt_vma_src(m))
        lndminjovrdmdni = log(dmt_grd(n  )/dmt_vma_src(m))
        ovr_src_snk_frc = 0.5_rkx * (erf(lndmaxjovrdmdni/sqrt2lngsdi) - &
                                     erf(lndminjovrdmdni/sqrt2lngsdi))
        ovr_src_snk_mss(m,n) = ovr_src_snk_frc * mss_frc_src(m)
      end do
    end do

    ! The following code from subroutine wnd_frc_thr_slt_get was placed
    ! here because tmp1 needs to be defined just once

    ryn_nbr_frc_thr_prx_opt = 0.38_rkx + 1331.0_rkx * (100.0_rkx*dmt_slt_opt)**1.56_rkx

    if (ryn_nbr_frc_thr_prx_opt < 0.03_rkx) then
      write(stderr,*) 'dstmbl: ryn_nbr_frc_thr_prx_opt < 0.03'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    else if (ryn_nbr_frc_thr_prx_opt < 10.0_rkx) then
      ryn_nbr_frc_thr_opt_fnc = -1.0_rkx + 1.928_rkx * &
              (ryn_nbr_frc_thr_prx_opt**0.0922_rkx)
      ryn_nbr_frc_thr_opt_fnc = 0.1291_rkx * 0.1291_rkx / ryn_nbr_frc_thr_opt_fnc
    else
      ryn_nbr_frc_thr_opt_fnc = 1.0_rkx - 0.0858_rkx * &
              exp(-0.0617_rkx*(ryn_nbr_frc_thr_prx_opt-10.0_rkx))
      ryn_nbr_frc_thr_opt_fnc = 0.120_rkx * 0.120_rkx * &
              ryn_nbr_frc_thr_opt_fnc * ryn_nbr_frc_thr_opt_fnc
    end if

    icf_fct = 1.0_rkx + 6.0e-7_rkx / (dns_slt * grav * (dmt_slt_opt**2.5_rkx))
    dns_fct = dns_slt * grav * dmt_slt_opt
    tmp1 = sqrt(icf_fct * dns_fct * ryn_nbr_frc_thr_opt_fnc)

    ! Set basin factor to 1 for now

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    do c = begc, endc
      l = clm3%g%l%c%landunit(c)
      if (.not. clm3%g%l%lakpoi(l)) then
        mbl_bsn_fct(c) = 1.0_rkx
      end if
    end do

    ! Introducing particle diameter. Needed by atm model and by dry dep model.
    ! Taken from Charlie Zender's subroutines dst_psd_ini, dst_sz_rsl,
    ! grd_mk (dstpsd.F90) and subroutine lgn_evl (psdlgn.F90)

    ! Charlie allows logarithmic or linear option for size distribution
    ! however, he hardwires the distribution to logarithmic in his code
    ! therefore, I take his logarithmic code only
    ! furthermore, if dst_nbr == 4, he overrides the automatic grid calculation
    ! he currently works with dst_nbr = 4, so I only take the relevant code
    ! if ndst ever becomes different from 4, must add call grd_mk (dstpsd.F90)
    ! as done in subroutine dst_psd_ini
    ! note that here ndst = dst_nbr

    ! Override automatic grid with preset grid if available

    if (ndst == 4) then
      do n = 1, ndst
        dmt_min(n) = dmt_grd(n)                       ![m] Max diameter in bin
        dmt_max(n) = dmt_grd(n+1)                     ![m] Min diameter in bin
        dmt_ctr(n) = 0.5_rkx * (dmt_min(n)+dmt_max(n)) ![m] Diameter at bin ctr
        dmt_dlt(n) = dmt_max(n)-dmt_min(n)            ![m] Width of size bin
      end do
    else
      write(stderr,*) 'Dustini error: ndst must equal to 4 with current code'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if !end if ndst == 4

    ! Bin physical properties

    gsd_anl = 2.0_rkx    ! [frc] Geometric std dev PaG77 p. 2080 Table1
    ln_gsd = log(gsd_anl)
    dns_aer = 2.5e3_rkx   ! [kg m-3] Aerosol density

    ! Set a fundamental statistic for each bin

    dmt_vma = 3.5000e-6_rkx ! [m] Mass median diameter analytic She84 p.75 Table1

    ! Compute analytic size statistics
    ! Convert mass median diameter to number median diameter (call vma2nma)

    dmt_nma = dmt_vma * exp(-3.0_rkx*ln_gsd*ln_gsd) ! [m]

    ! Compute resolved size statistics for each size distribution
    ! In C. Zender's code call dst_sz_rsl

    do n = 1, ndst

      series_ratio = (dmt_max(n)/dmt_min(n))**(1.0_rkx/sz_nbr)
      sz_min(1) = dmt_min(n)
      do m = 2, sz_nbr                            ! Loop starts at 2
        sz_min(m) = sz_min(m-1) * series_ratio
      end do

      ! Derived grid values
      do m = 1, sz_nbr-1                          ! Loop ends at sz_nbr-1
        sz_max(m) = sz_min(m+1)                  ! [m]
      end do
      sz_max(sz_nbr) = dmt_max(n)                 ! [m]

      ! Final derived grid values
      do m = 1, sz_nbr
        sz_ctr(m) = 0.5_rkx * (sz_min(m)+sz_max(m))
        sz_dlt(m) = sz_max(m)-sz_min(m)
      end do

      lngsdsqrttwopi_rcp = 1.0_rkx / (ln_gsd*sqrt(2.0_rkx*rpi))
      dmt_vwr(n) = 0.0_rkx ! [m] Mass wgted diameter resolved
      vlm_rsl(n) = 0.0_rkx ! [m3 m-3] Volume concentration resolved

      do m = 1, sz_nbr

        ! Evaluate lognormal distribution for these sizes (call lgn_evl)
        tmp = log(sz_ctr(m)/dmt_nma) / ln_gsd
        lgn_dst = lngsdsqrttwopi_rcp * exp(-0.5_rkx*tmp*tmp) / sz_ctr(m)

        ! Integrate moments of size distribution
        dmt_vwr(n) = dmt_vwr(n) + sz_ctr(m) *                    &
               rpi / 6.0_rkx * (sz_ctr(m)**3.0_rkx) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn
        vlm_rsl(n) = vlm_rsl(n) +                                &
               rpi / 6.0_rkx * (sz_ctr(m)**3.0_rkx) * & ![m3] Volume
               lgn_dst * sz_dlt(m)                ![# m-3] Number concentrn

      end do
      dmt_vwr(n) = dmt_vwr(n) / vlm_rsl(n) ![m] Mass weighted diameter resolved
    end do

    ! calculate correction to Stokes' settling velocity (subroutine stk_crc_get)

    eps_max = 1.0e-4_rkx
    dns_mdp = 100000._rkx / (295.0_rkx*rair) ![kg m-3] const prs_mdp & tpt_vrt

    ! Size-independent thermokinetic properties

    vsc_dyn_atm = 1.72e-5_rkx * ((295.0_rkx/273.0_rkx)**1.5_rkx) * 393.0_rkx / &
         (295.0_rkx+120.0_rkx)      ![kg m-1 s-1] RoY94 p.102 tpt_mdp=295.0
    mfp_atm = 2.0_rkx * vsc_dyn_atm / &  !SeP97 p. 455 constant prs_mdp, tpt_mdp
         (100000._rkx*sqrt(8.0_rkx/(rpi*rair*295.0_rkx)))
    vsc_knm_atm = vsc_dyn_atm / dns_mdp ![m2 s-1] Kinematic viscosity of air

    do m = 1, ndst
      slp_crc(m) = 1.0_rkx + 2.0_rkx * mfp_atm *                      &
          (1.257_rkx+0.4_rkx*exp(-1.1_rkx*dmt_vwr(m)/(2.0_rkx*mfp_atm))) / &
          dmt_vwr(m)  ! [frc] Slip correction factor SeP97 p.464
      vlc_stk(m) = (1.0_rkx/18.0_rkx) * dmt_vwr(m) * dmt_vwr(m) * dns_aer * &
          grav * slp_crc(m) / vsc_dyn_atm ! [m s-1] SeP97 p.466
    end do

    ! For Reynolds number flows Re < 0.1 Stokes' velocity is valid for
    ! vlc_grv SeP97 p. 466 (8.42). For larger Re, inertial effects become
    ! important and empirical drag coefficients must be employed
    ! Implicit equation for Re, Cd, and Vt is SeP97 p. 467 (8.44)
    ! Using Stokes' velocity rather than iterative solution with empirical
    ! drag coefficient causes 60% errors for D = 200 um SeP97 p. 468

    ! Iterative solution for drag coefficient, Reynolds number,
    ! and terminal veloc
    do m = 1, ndst

      ! Initialize accuracy and counter
      eps_crr = eps_max + 1.0_rkx  ![frc] Current relative accuracy
      itr_idx = 0                 ![idx] Counting index

      ! Initial guess for vlc_grv is exact for Re < 0.1
      vlc_grv(m) = vlc_stk(m)     ![m s-1]

      do while(eps_crr > eps_max)

        ! Save terminal velocity for convergence test
        vlc_grv_old = vlc_grv(m) ![m s-1]
        ryn_nbr_grv(m) = vlc_grv(m) * dmt_vwr(m) / vsc_knm_atm !SeP97 p.460

        ! Update drag coefficient based on new Reynolds number
        if (ryn_nbr_grv(m) < 0.1_rkx) then
          !Stokes' law Sep97 p.463 (8.32)
          cff_drg_grv(m) = 24.0_rkx / ryn_nbr_grv(m)
        else if (ryn_nbr_grv(m) < 2.0_rkx) then
          cff_drg_grv(m) = (24.0_rkx/ryn_nbr_grv(m)) *    &
                 (1.0_rkx + 3.0_rkx*ryn_nbr_grv(m)/16.0_rkx + &
               9.0_rkx*ryn_nbr_grv(m)*ryn_nbr_grv(m)*     &
                log(2.0_rkx*ryn_nbr_grv(m))/160.0_rkx)        !Sep97 p.463 (8.32)
        else if (ryn_nbr_grv(m) < 500.0_rkx) then
          cff_drg_grv(m) = (24.0_rkx/ryn_nbr_grv(m)) * &
               (1.0_rkx + 0.15_rkx*ryn_nbr_grv(m)**0.687_rkx) !Sep97 p.463 (8.32)
        else if (ryn_nbr_grv(m) < 2.0e5_rkx) then
          cff_drg_grv(m) = 0.44_rkx                         !Sep97 p.463 (8.32)
        else
          write(stderr,'(a,es9.2)') "ryn_nbr_grv(m) = ",ryn_nbr_grv(m)
          write(stderr,*) &
                  'Dustini error: Reynolds number too large in stk_crc_get()'
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if

        ! Update terminal velocity based on new Reynolds number and drag coeff
        ! [m s-1] Terminal veloc SeP97 p.467 (8.44)

        vlc_grv(m) = sqrt(4.0_rkx * grav * dmt_vwr(m) * slp_crc(m) * dns_aer / &
             (3.0_rkx*cff_drg_grv(m)*dns_mdp))
        eps_crr = abs((vlc_grv(m)-vlc_grv_old)/vlc_grv(m)) !Relative convergence
        if (itr_idx == 12) then
          ! Numerical pingpong may occur when Re = 0.1, 2.0, or 500.0
          ! due to discontinuities in derivative of drag coefficient
          vlc_grv(m) = 0.5_rkx * (vlc_grv(m)+vlc_grv_old)  ! [m s-1]
        end if
        if (itr_idx > 20) then
          write(stderr,*) 'Dustini error: Terminal velocity not converging ',&
                  ' in stk_crc_get(), breaking loop...'
          exit
        end if
        itr_idx = itr_idx + 1

      end do
    end do   !end loop over size

    ! Compute factors to convert Stokes' settling velocities to
    ! actual settling velocities

    do m = 1, ndst
      stk_crc(m) = vlc_grv(m) / vlc_stk(m)
    end do

  end subroutine Dustini

end module mod_clm_dust
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
