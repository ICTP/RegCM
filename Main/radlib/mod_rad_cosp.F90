!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!
! Satellite simulators of COSP2 (Bodas-Salcedo et al., 2011; Swales et al.,
! 2018, external/cosp2) on the cloud seen by the radiation (icosp = 1):
! what the ISCCP and MODIS products would be for the model state.
!
! At every radiation call the layer cloud fraction, in cloud water path,
! ice fraction and effective radii of the radiation are split in ncol
! subcolumns (SCOPS, maximum random overlap, as the clh, clm and cll
! diagnostic), with the 0.67 micron optical depth of the CCM visible band
! (liquid Slingo, ice Ebert and Curry), the 11 micron emissivity of the CCM
! long wave (Savijarvi and Raisanen) and the MODIS near infrared single
! scattering properties of COSP. The simulator results are averaged over
! the output interval: the cloud fractions and brightness temperatures over
! the valid calls (the passive sensor products are day time only), the
! cloud properties weighted by the cloud fraction they refer to, as the
! monthly means of the satellite products.
!
module mod_rad_cosp

  use mod_intkinds, only : ik4, ik8
  use mod_realkinds, only : rk8, rkx
  use mod_dynparam, only : kz, kzp1, jci1, jci2, ici1, ici2
  use mod_runparams, only : ipptls, iqv, iqc, iqi
  use mod_regcm_types, only : mod_2_rad
  use mod_rad_common, only : cldwpmin, cloud_fice, cloud_rel_cam, &
    cloud_rei
  use mod_constants, only : egrav, rgas
  use mod_message, only : fatal
  use mod_memutil, only : getmem
  use mod_mppparam, only : global_dot_jstart, global_dot_istart
  use mod_runparams, only : rcmtimer
  use mod_outvars, only : ncospvars, rad_cosp_out
  use mod_rad_radiation, only : abarl, bbarl, abari, bbari
  use cosp_kinds, only : wp
  use mod_cosp_config, only : parasol_nrefl
  use mod_cosp_stats, only : cosp_optical_inputs, cosp_column_inputs, &
    radar_cfg
  use mod_cosp_rttov_util, only : rttov_cfg
  use mod_cosp, only : cosp_init, cosp_simulator, cosp_outputs
  use mod_scops, only : scops
  use mod_rng, only : rng_state, init_rng
  use cosp_optics, only : cosp_simulator_optics, modis_optics

  implicit none

  private

  public :: allocate_mod_rad_cosp, cosp_driver, cosp_output

  ! Number of subcolumns
  integer(ik4), parameter :: ncol = 20
  ! Overlap of the cloud layers in SCOPS: 3 maximum random
  integer(ik4), parameter :: overlap = 3
  ! Columns of a call to the simulators
  integer(ik4), parameter :: nchunk = 128
  ! Missing values of the simulators (R_UNDEF, ICARUS output_missing_value
  ! -1.0e30 in single or double precision) are below undef
  real(wp), parameter :: undef = -1.0e29_wp
  ! Surface emissivity at 11 micron
  real(wp), parameter :: emsfc_lw = 0.99_wp

  ! Field of the cloud fraction weighting each output field (0: none)
  integer(ik4), dimension(ncospvars), parameter :: iweight = &
    [ 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 7, 8, 9, 8, 9, 7, 8, 9 ]

  ! Sums of the field (times the weight) and of the weight (or counts)
  real(rk8), pointer, contiguous, dimension(:,:,:) :: csum => null( )
  real(rk8), pointer, contiguous, dimension(:,:,:) :: cwgt => null( )

  type(cosp_outputs), save :: cospout
  integer(ik4) :: npr

  contains

  subroutine allocate_mod_rad_cosp
    implicit none
    type(radar_cfg) :: rcfg
    type(rttov_cfg), dimension(:), allocatable :: rttov_configs
    character(len=256), dimension(0) :: rttov_namelists
    character(len=64) :: micro_scheme
    logical :: lrttov
    npr = (jci2-jci1+1)*(ici2-ici1+1)
    call getmem(csum,jci1,jci2,ici1,ici2,1,ncospvars,'cosp:csum')
    call getmem(cwgt,jci1,jci2,ici1,ici2,1,ncospvars,'cosp:cwgt')
    csum(:,:,:) = 0.0_rk8
    cwgt(:,:,:) = 0.0_rk8
    lrttov = .false.
    micro_scheme = 'MMF_v3_single_moment'
    ! ISCCP and MODIS only; ISCCP cloud top from the 10.5 micron
    ! brightness temperature (1), lowest level of that temperature (2)
    call cosp_init(.true.,.true.,.false.,.false.,.false.,.false.,  &
                   .false.,.false.,lrttov,94.0_wp,-1.0_wp,0,0,1,2,0, &
                   rcfg,.false.,.false.,0,kz,micro_scheme,0,         &
                   rttov_namelists,rttov_configs)
    ! ISCCP subcolumn cloud top pressure and optical depth, used by MODIS
    allocate(cospout%isccp_boxtau(npr,ncol))
    allocate(cospout%isccp_boxptop(npr,ncol))
    allocate(cospout%isccp_totalcldarea(npr))
    allocate(cospout%isccp_meanptop(npr))
    allocate(cospout%isccp_meantaucld(npr))
    allocate(cospout%isccp_meanalbedocld(npr))
    allocate(cospout%isccp_meantb(npr))
    allocate(cospout%isccp_meantbclr(npr))
    allocate(cospout%modis_cloud_fraction_total_mean(npr))
    allocate(cospout%modis_cloud_fraction_water_mean(npr))
    allocate(cospout%modis_cloud_fraction_ice_mean(npr))
    allocate(cospout%modis_cloud_fraction_high_mean(npr))
    allocate(cospout%modis_cloud_fraction_mid_mean(npr))
    allocate(cospout%modis_cloud_fraction_low_mean(npr))
    allocate(cospout%modis_optical_thickness_total_mean(npr))
    allocate(cospout%modis_optical_thickness_water_mean(npr))
    allocate(cospout%modis_optical_thickness_ice_mean(npr))
    allocate(cospout%modis_cloud_particle_size_water_mean(npr))
    allocate(cospout%modis_cloud_particle_size_ice_mean(npr))
    allocate(cospout%modis_cloud_top_pressure_total_mean(npr))
    allocate(cospout%modis_liquid_water_path_mean(npr))
    allocate(cospout%modis_ice_water_path_mean(npr))
  end subroutine allocate_mod_rad_cosp
  !
  ! Run the simulators on the cloud seen by the radiation, the same for
  ! both radiation schemes: cloud fraction and in cloud water of the
  ! layers with a water path of at least cldwpmin, ice fraction and
  ! effective radii of mod_rad_common (without the sulfate indirect
  ! effect on the liquid radius); add the results to the sums.
  !
  subroutine cosp_driver(m2r)
    implicit none
    type(mod_2_rad), intent(inout) :: m2r
    real(rk8) :: tk, ql, qi, cf, wp_in, fi, re_l, re_i, dzk, iwc
    type(cosp_optical_inputs) :: cin
    type(cosp_column_inputs) :: cst
    type(rng_state), dimension(:), allocatable :: rngs
    integer(ik4), dimension(:), allocatable :: seed
    real(wp), dimension(:,:), allocatable :: tca, conv, tl, ti, dem, &
      rl, ri, rho
    real(wp), dimension(:,:,:), allocatable :: taul, taui, sizel, sizei
    character(len=256), dimension(100) :: status
    integer(ik4) :: n0, n1, n2, np, p, n, k, nj, j, i
    integer(ik8) :: iseed
    real(rk8) :: zh, zb, kabs, kabsl, kabsi
    nj = jci2-jci1+1
    ! The simulators run on the host: copy back the inputs from the
    ! device (nothing to do with unified memory)
!$acc update self(m2r%tatms,m2r%phatms,m2r%pfatms,m2r%qxatms) if_present
!$acc update self(m2r%cldfrc,m2r%cldlwc,m2r%deltaz,m2r%coszrs) if_present
!$acc update self(m2r%tg,m2r%ldmsk,m2r%ht,m2r%xlat) if_present
    do n0 = 1, npr, nchunk
      n1 = n0
      n2 = min(npr,n0+nchunk-1)
      np = n2-n1+1
      call construct_inputs(np,cin,cst)
      allocate(rngs(np),seed(np))
      allocate(tca(np,kz),conv(np,kz),tl(np,kz),ti(np,kz),dem(np,kz))
      allocate(rl(np,kz),ri(np,kz),rho(np,kz))
      allocate(taul(np,ncol,kz),taui(np,ncol,kz))
      allocate(sizel(np,ncol,kz),sizei(np,ncol,kz))
      do p = 1, np
        n = n1+p-1
        j = jci1 + mod(n-1,nj)
        i = ici1 + (n-1)/nj
        cst%sunlit(p) = merge(1,0,m2r%coszrs(j,i) >= 1.0e-3_rkx)
        cst%skt(p) = m2r%tg(j,i)
        cst%land(p) = merge(1.0_wp,0.0_wp,m2r%ldmsk(j,i) == 1)
        cst%rttov_sfcmask(p) = cst%land(p)
        cst%surfelev(p) = m2r%ht(j,i)/egrav
        cst%psfc(p) = m2r%pfatms(j,i,kzp1)
        cst%lat(p) = 0.0_wp
        cst%lon(p) = 0.0_wp
        cst%u_sfc(p) = 0.0_wp
        cst%v_sfc(p) = 0.0_wp
        cst%sza(p) = 0.0_wp
        cst%rttov_date(p,:) = 0.0_wp
        cst%rttov_time(p,:) = 0.0_wp
        do k = 1, kzp1
          cst%phalf(p,k) = m2r%pfatms(j,i,k)
        end do
        ! Heights above the sea level (km): layer middle and bottom
        zb = cst%surfelev(p)
        do k = kz, 1, -1
          zh = zb + m2r%deltaz(j,i,k)
          cst%hgt_matrix_half(p,k) = zb*0.001_wp
          cst%hgt_matrix(p,k) = 0.5_wp*(zb+zh)*0.001_wp
          zb = zh
        end do
        do k = 1, kz
          tk = real(m2r%tatms(j,i,k),rk8)
          dzk = real(m2r%deltaz(j,i,k),rk8)
          cst%at(p,k) = tk
          cst%pfull(p,k) = m2r%phatms(j,i,k)
          ! Specific humidity from the mixing ratio
          cst%qv(p,k) = m2r%qxatms(j,i,k,iqv) / &
                        (1.0_rk8+m2r%qxatms(j,i,k,iqv))
          rho(p,k) = m2r%phatms(j,i,k)/(rgas*tk)
          ! Cloud of the radiation (colmod3, rrtmg_driver)
          wp_in = real(m2r%cldlwc(j,i,k),rk8)*dzk
          cf = 0.0_rk8
          if ( wp_in > cldwpmin ) then
            cf = min(max(real(m2r%cldfrc(j,i,k),rk8),0.0_rk8),1.0_rk8)
          end if
          ql = real(m2r%qxatms(j,i,k,iqc),rk8)
          if ( ipptls > 1 ) then
            qi = real(m2r%qxatms(j,i,k,iqi),rk8)
            fi = cloud_fice(.true.,ql,qi,tk)
          else
            fi = cloud_fice(.false.,ql,0.0_rk8,tk)
          end if
          iwc = wp_in*fi/dzk
          re_l = cloud_rel_cam(tk,m2r%ldmsk(j,i) == 1)
          re_i = cloud_rei(iwc,tk,real(m2r%xlat(j,i),rk8))
          if ( cf > 0.0_rk8 ) then
            tca(p,k) = cf
            ! In cloud 0.67 micron optical depth, CCM visible band
            tl(p,k) = wp_in*(1.0_rk8-fi)*(abarl(1)+bbarl(1)/re_l)
            ti(p,k) = wp_in*fi*(abari(1)+bbari(1)/re_i)
            ! In cloud 11 micron emissivity, as the CCM long wave
            kabsi = 0.04_rk8 + 0.08_rk8*exp(-0.018_rk8*re_i)
            kabsl = 0.31_rk8*exp(-0.08_rk8*re_l)
            kabs = kabsl*(1.0_rk8-fi) + kabsi*fi
            dem(p,k) = 1.0_rk8 - exp(-min(kabs*wp_in,25.0_rk8))
          else
            tca(p,k) = 0.0_wp
            tl(p,k) = 0.0_wp
            ti(p,k) = 0.0_wp
            dem(p,k) = 0.0_wp
          end if
          rl(p,k) = re_l
          ri(p,k) = re_i
          ! Grid mean cloud water and ice (kg/kg)
          cst%cloudliq(p,k) = wp_in*tca(p,k)*(1.0_rk8-fi) / &
                              (rho(p,k)*dzk)*0.001_wp
          cst%cloudice(p,k) = wp_in*tca(p,k)*fi / &
                              (rho(p,k)*dzk)*0.001_wp
          cst%deffliq(p,k) = 2.0_wp*re_l
          cst%deffice(p,k) = 2.0_wp*re_i
        end do
        cst%tca(p,:) = tca(p,:)
        ! Seed of the subcolumns: global point and time step, the same
        ! for any domain decomposition
        iseed = int(global_dot_jstart+j-1,ik8)*7919_ik8 + &
                int(global_dot_istart+i-1,ik8)*104729_ik8 + &
                rcmtimer%lcount*15485863_ik8
        seed(p) = int(1_ik8+mod(iseed,2147483646_ik8),ik4)
      end do
      cst%o3(:,:) = 0.0_wp
      cst%co2(:,:) = 5.241e-04_wp
      cst%ch4(:,:) = 9.139e-07_wp
      cst%n2o(:,:) = 4.665e-07_wp
      cst%co(:,:) = 2.098e-07_wp
      cst%so2(:,:) = 2.0e-11_wp
      cst%fl_rain(:,:) = 0.0_wp
      cst%fl_snow(:,:) = 0.0_wp
      conv(:,:) = 0.0_wp
      call init_rng(rngs,seed)
      call scops(np,kz,ncol,rngs,tca,conv,overlap,cin%frac_out,0)
      call cosp_simulator_optics(np,ncol,kz,cin%frac_out,tl+ti,tl+ti, &
                                 cin%tau_067)
      call cosp_simulator_optics(np,ncol,kz,cin%frac_out,dem,dem, &
                                 cin%emiss_11)
      call cosp_simulator_optics(np,ncol,kz,cin%frac_out,tl,tl,taul)
      call cosp_simulator_optics(np,ncol,kz,cin%frac_out,ti,ti,taui)
      call cosp_simulator_optics(np,ncol,kz,cin%frac_out,rl,rl,sizel)
      call cosp_simulator_optics(np,ncol,kz,cin%frac_out,ri,ri,sizei)
      call modis_optics(np,kz,ncol,taul,sizel,taui,sizei, &
                        cin%fracliq,cin%asym,cin%ss_alb)
      cin%emsfc_lw = emsfc_lw
      status = cosp_simulator(cin,cst,cospout,n1,n2,.false.)
      if ( any(len_trim(status) > 0) ) then
        call fatal(__FILE__,__LINE__,trim(status(1)))
      end if
      call destroy_inputs(cin,cst)
      deallocate(rngs,seed,tca,conv,tl,ti,dem,rl,ri,rho)
      deallocate(taul,taui,sizel,sizei)
    end do
    call accumulate( 1,cospout%isccp_totalcldarea)
    call accumulate( 2,cospout%isccp_meanptop)
    call accumulate( 3,cospout%isccp_meantaucld)
    call accumulate( 4,cospout%isccp_meanalbedocld)
    call accumulate( 5,cospout%isccp_meantb)
    call accumulate( 6,cospout%isccp_meantbclr)
    call accumulate( 7,cospout%modis_cloud_fraction_total_mean)
    call accumulate( 8,cospout%modis_cloud_fraction_water_mean)
    call accumulate( 9,cospout%modis_cloud_fraction_ice_mean)
    call accumulate(10,cospout%modis_cloud_fraction_high_mean)
    call accumulate(11,cospout%modis_cloud_fraction_mid_mean)
    call accumulate(12,cospout%modis_cloud_fraction_low_mean)
    call accumulate(13,cospout%modis_optical_thickness_total_mean)
    call accumulate(14,cospout%modis_optical_thickness_water_mean)
    call accumulate(15,cospout%modis_optical_thickness_ice_mean)
    call accumulate(16,cospout%modis_cloud_particle_size_water_mean)
    call accumulate(17,cospout%modis_cloud_particle_size_ice_mean)
    call accumulate(18,cospout%modis_cloud_top_pressure_total_mean)
    call accumulate(19,cospout%modis_liquid_water_path_mean)
    call accumulate(20,cospout%modis_ice_water_path_mean)

    contains

    ! Add a field, valid where it is not undefined, weighted by the
    ! cloud fraction field iweight (in %) when it has one.
    subroutine accumulate(iv,f)
      implicit none
      integer(ik4), intent(in) :: iv
      real(wp), dimension(:), intent(in) :: f
      integer(ik4) :: m
      real(rk8) :: w
      real(wp), pointer, dimension(:) :: fw
      do n = 1, npr
        if ( f(n) < undef ) cycle
        ! MODIS gives a zero cloud top pressure to cloudy subcolumns
        ! without a retrieved cloud top
        if ( (iv == 2 .or. iv == 18) .and. f(n) <= 0.0_wp ) cycle
        j = jci1 + mod(n-1,nj)
        i = ici1 + (n-1)/nj
        if ( iweight(iv) > 0 ) then
          m = iweight(iv)
          fw => weight_field(m)
          if ( fw(n) < undef ) cycle
          w = fw(n)*0.01_rk8
        else
          w = 1.0_rk8
        end if
        csum(j,i,iv) = csum(j,i,iv) + f(n)*w
        cwgt(j,i,iv) = cwgt(j,i,iv) + w
      end do
    end subroutine accumulate

  end subroutine cosp_driver
  !
  ! Cloud fraction output field m of the weighted means
  !
  function weight_field(m) result(fw)
    implicit none
    integer(ik4), intent(in) :: m
    real(wp), pointer, dimension(:) :: fw
    select case (m)
      case (1)
        fw => cospout%isccp_totalcldarea
      case (7)
        fw => cospout%modis_cloud_fraction_total_mean
      case (8)
        fw => cospout%modis_cloud_fraction_water_mean
      case (9)
        fw => cospout%modis_cloud_fraction_ice_mean
      case default
        call fatal(__FILE__,__LINE__,'No weight field')
    end select
  end function weight_field
  !
  ! Means over the output interval to the output fields, missing value
  ! where there is no valid value; restart the sums.
  !
  subroutine cosp_output(rmiss)
    implicit none
    real(rkx), intent(in) :: rmiss
    integer(ik4) :: iv, i, j
    do iv = 1, ncospvars
      if ( .not. associated(rad_cosp_out(iv)%p) ) cycle
      do i = ici1, ici2
        do j = jci1, jci2
          if ( cwgt(j,i,iv) > 0.0_rk8 ) then
            rad_cosp_out(iv)%p(j,i) = real(csum(j,i,iv)/cwgt(j,i,iv),rkx)
          else
            rad_cosp_out(iv)%p(j,i) = rmiss
          end if
        end do
      end do
    end do
    csum(:,:,:) = 0.0_rk8
    cwgt(:,:,:) = 0.0_rk8
  end subroutine cosp_output

  subroutine construct_inputs(np,cin,cst)
    implicit none
    integer(ik4), intent(in) :: np
    type(cosp_optical_inputs), intent(out) :: cin
    type(cosp_column_inputs), intent(out) :: cst
    cin%npoints = np
    cin%ncolumns = ncol
    cin%nlevels = kz
    cin%ninst_rttov = 0
    cin%npart = 4
    cin%nrefl = parasol_nrefl
    allocate(cin%frac_out(np,ncol,kz))
    allocate(cin%tau_067(np,ncol,kz))
    allocate(cin%emiss_11(np,ncol,kz))
    allocate(cin%fracliq(np,ncol,kz))
    allocate(cin%asym(np,ncol,kz))
    allocate(cin%ss_alb(np,ncol,kz))
    cst%npoints = np
    cst%ncolumns = ncol
    cst%nlevels = kz
    allocate(cst%sunlit(np),cst%skt(np),cst%land(np),cst%psfc(np))
    allocate(cst%at(np,kz),cst%pfull(np,kz),cst%phalf(np,kzp1))
    allocate(cst%qv(np,kz),cst%o3(np,kz),cst%hgt_matrix(np,kz))
    allocate(cst%hgt_matrix_half(np,kz))
    allocate(cst%u_sfc(np),cst%v_sfc(np),cst%lat(np),cst%lon(np))
    allocate(cst%rttov_sfcmask(np),cst%surfelev(np),cst%sza(np))
    allocate(cst%co(np,kz),cst%n2o(np,kz),cst%ch4(np,kz))
    allocate(cst%co2(np,kz),cst%so2(np,kz))
    allocate(cst%cloudice(np,kz),cst%cloudliq(np,kz))
    allocate(cst%deffliq(np,kz),cst%deffice(np,kz))
    allocate(cst%fl_snow(np,kz),cst%fl_rain(np,kz),cst%tca(np,kz))
    allocate(cst%rttov_date(np,3),cst%rttov_time(np,3))
  end subroutine construct_inputs

  subroutine destroy_inputs(cin,cst)
    implicit none
    type(cosp_optical_inputs), intent(inout) :: cin
    type(cosp_column_inputs), intent(inout) :: cst
    deallocate(cin%frac_out,cin%tau_067,cin%emiss_11)
    deallocate(cin%fracliq,cin%asym,cin%ss_alb)
    deallocate(cst%sunlit,cst%skt,cst%land,cst%psfc)
    deallocate(cst%at,cst%pfull,cst%phalf,cst%qv,cst%o3)
    deallocate(cst%hgt_matrix,cst%hgt_matrix_half)
    deallocate(cst%u_sfc,cst%v_sfc,cst%lat,cst%lon)
    deallocate(cst%rttov_sfcmask,cst%surfelev,cst%sza)
    deallocate(cst%co,cst%n2o,cst%ch4,cst%co2,cst%so2)
    deallocate(cst%cloudice,cst%cloudliq,cst%deffliq,cst%deffice)
    deallocate(cst%fl_snow,cst%fl_rain,cst%tca)
    deallocate(cst%rttov_date,cst%rttov_time)
  end subroutine destroy_inputs

end module mod_rad_cosp
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
