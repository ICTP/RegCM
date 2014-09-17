!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
module mod_ncout
!
  use mod_intkinds
  use mod_realkinds
  use mod_memutil
  use mod_dynparam
  use mod_runparams
  use mod_mppparam
  use mod_ncstream_types
  use mod_ncstream
  use mod_outvars
  use mod_service
  use netcdf

  implicit none

  public :: init_output_streams
  public :: dispose_output_streams
  public :: write_record_output_stream
  public :: writevar_output_stream

  type varspan
    integer(ik4) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
  end type varspan

  type regcm_stream
    type(nc_output_stream) , pointer , dimension(:) :: ncout
    character(len=8) , pointer , dimension(:) :: cname_base
    type(ncoutstream_params) :: opar
    logical :: l_sub = .false.
    integer(ik4) :: nvar = 0
    integer(ik4) :: nfiles = 0
    type(nc_varlist) :: ncvars
    integer(ik4) :: jg1 , jg2 , ig1 , ig2
    integer(ik4) :: jl1 , jl2 , il1 , il2
  end type regcm_stream

  integer(ik4) , parameter :: nbase = 5

  integer(ik4) , parameter :: natm2dvars = 4 + nbase
  integer(ik4) , parameter :: natm3dvars = 55
  integer(ik4) , parameter :: natmvars = natm2dvars+natm3dvars

  integer(ik4) , parameter :: nsrf2dvars = 18 + nbase
  integer(ik4) , parameter :: nsrf3dvars = 7
  integer(ik4) , parameter :: nsrfvars = nsrf2dvars+nsrf3dvars

  integer(ik4) , parameter :: nsts2dvars = 7 + nbase
  integer(ik4) , parameter :: nsts3dvars = 5
  integer(ik4) , parameter :: nstsvars = nsts2dvars+nsts3dvars

  integer(ik4) , parameter :: nsub2dvars = 6 + nbase
  integer(ik4) , parameter :: nsub3dvars = 6
  integer(ik4) , parameter :: nsubvars = nsub2dvars+nsub3dvars

  integer(ik4) , parameter :: nlak2dvars = 12 + nbase
  integer(ik4) , parameter :: nlak3dvars = 1
  integer(ik4) , parameter :: nlakvars = nlak2dvars+nlak3dvars

  integer(ik4) , parameter :: nrad2dvars = 12 + nbase
  integer(ik4) , parameter :: nrad3dvars = 4
  integer(ik4) , parameter :: nrad4dvars = 2
  integer(ik4) , parameter :: nradvars = nrad2dvars+nrad3dvars+nrad4dvars

  integer(ik4) , parameter :: nopt2dvars = 5 + nbase
  integer(ik4) , parameter :: nopt3dvars = 3
  integer(ik4) , parameter :: noptvars = nopt2dvars+nopt3dvars

  integer(ik4) , parameter :: nche2dvars = 7 + nbase
  integer(ik4) , parameter :: nche3dvars = 12
  integer(ik4) , parameter :: nchevars = nche2dvars+nche3dvars

  integer(ik4) , parameter :: nslaboc2dvars = nbase
  integer(ik4) , parameter :: nslaboc3dvars = 1
  integer(ik4) , parameter :: nslabocvars = nslaboc2dvars + nslaboc3dvars

  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_atm => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_atm => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_srf => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_srf => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_sts => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_sts => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_sub => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_sub => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_lak => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_lak => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_rad => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_rad => null()
  type(ncvariable4d_real) , save , pointer , &
    dimension(:) :: v4dvar_rad => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_opt => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_opt => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_che => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_che => null()
  type(ncvariable2d_real) , save , pointer , &
    dimension(:) :: v2dvar_slaboc => null()
  type(ncvariable3d_real) , save , pointer , &
    dimension(:) :: v3dvar_slaboc => null()

  integer(ik4) :: maxstreams

  logical :: parallel_out

  integer(ik4) , public :: atm_stream = -1
  integer(ik4) , public :: srf_stream = -1
  integer(ik4) , public :: sub_stream = -1
  integer(ik4) , public :: rad_stream = -1
  integer(ik4) , public :: lak_stream = -1
  integer(ik4) , public :: sts_stream = -1
  integer(ik4) , public :: opt_stream = -1
  integer(ik4) , public :: che_stream = -1
  integer(ik4) , public :: slaboc_stream = -1
!
  type(regcm_stream) , pointer , save , dimension(:) :: outstream

  logical , public , dimension(natmvars) :: enable_atm_vars
  logical , public , dimension(nsrfvars) :: enable_srf_vars
  logical , public , dimension(nstsvars) :: enable_sts_vars
  logical , public , dimension(nsubvars) :: enable_sub_vars
  logical , public , dimension(nlakvars) :: enable_lak_vars
  logical , public , dimension(nradvars) :: enable_rad_vars
  logical , public , dimension(noptvars) :: enable_opt_vars
  logical , public , dimension(nchevars) :: enable_che_vars

  integer(ik4) , parameter :: atm_xlon  = 1
  integer(ik4) , parameter :: atm_xlat  = 2
  integer(ik4) , parameter :: atm_mask  = 3
  integer(ik4) , parameter :: atm_topo  = 4
  integer(ik4) , parameter :: atm_ps    = 5
  integer(ik4) , parameter :: atm_tpr   = 6
  integer(ik4) , parameter :: atm_tsn   = 7
  integer(ik4) , parameter :: atm_tgb   = 8
  integer(ik4) , parameter :: atm_tsw   = 9

  integer(ik4) , parameter :: atm_u            = 1
  integer(ik4) , parameter :: atm_v            = 2
  integer(ik4) , parameter :: atm_t            = 3
  integer(ik4) , parameter :: atm_omega        = 4
  integer(ik4) , parameter :: atm_qv           = 5
  integer(ik4) , parameter :: atm_qc           = 6
  integer(ik4) , parameter :: atm_rh           = 7
  integer(ik4) , parameter :: atm_qr           = 8
  integer(ik4) , parameter :: atm_qi           = 9
  integer(ik4) , parameter :: atm_qs           = 10
  integer(ik4) , parameter :: atm_zf           = 11
  integer(ik4) , parameter :: atm_zh           = 12
  integer(ik4) , parameter :: atm_pf           = 13
  integer(ik4) , parameter :: atm_ph           = 14
  integer(ik4) , parameter :: atm_q_detr       = 15
  integer(ik4) , parameter :: atm_rainls       = 16
  integer(ik4) , parameter :: atm_raincc       = 17
  integer(ik4) , parameter :: atm_tke          = 18
  integer(ik4) , parameter :: atm_kth          = 19
  integer(ik4) , parameter :: atm_kzm          = 20
  integer(ik4) , parameter :: atm_tten_adh     = 21
  integer(ik4) , parameter :: atm_tten_adv     = 22
  integer(ik4) , parameter :: atm_tten_tbl     = 23
  integer(ik4) , parameter :: atm_tten_dif     = 24
  integer(ik4) , parameter :: atm_tten_bdy     = 25
  integer(ik4) , parameter :: atm_tten_con     = 26
  integer(ik4) , parameter :: atm_tten_adi     = 27
  integer(ik4) , parameter :: atm_tten_rad     = 28
  integer(ik4) , parameter :: atm_tten_lsc     = 29
  integer(ik4) , parameter :: atm_qten_adh     = 30
  integer(ik4) , parameter :: atm_qten_adv     = 31
  integer(ik4) , parameter :: atm_qten_tbl     = 32
  integer(ik4) , parameter :: atm_qten_dif     = 33
  integer(ik4) , parameter :: atm_qten_bdy     = 34
  integer(ik4) , parameter :: atm_qten_con     = 35
  integer(ik4) , parameter :: atm_qten_adi     = 36
  integer(ik4) , parameter :: atm_qten_rad     = 37
  integer(ik4) , parameter :: atm_qten_lsc     = 38
  integer(ik4) , parameter :: atm_stats_supw   = 39
  integer(ik4) , parameter :: atm_stats_supc   = 40
  integer(ik4) , parameter :: atm_stats_detw   = 41
  integer(ik4) , parameter :: atm_stats_detc   = 42
  integer(ik4) , parameter :: atm_stats_erow   = 43
  integer(ik4) , parameter :: atm_stats_eroc   = 44
  integer(ik4) , parameter :: atm_stats_evw    = 45
  integer(ik4) , parameter :: atm_stats_evc    = 46
  integer(ik4) , parameter :: atm_stats_con1w  = 47
  integer(ik4) , parameter :: atm_stats_con1c  = 48
  integer(ik4) , parameter :: atm_stats_con2w  = 49
  integer(ik4) , parameter :: atm_stats_con2c  = 50
  integer(ik4) , parameter :: atm_stats_dep    = 51
  integer(ik4) , parameter :: atm_stats_melt   = 52
  integer(ik4) , parameter :: atm_stats_frz    = 53
  integer(ik4) , parameter :: atm_stats_rainev = 54
  integer(ik4) , parameter :: atm_stats_snowev = 55

  integer(ik4) , parameter :: srf_xlon   = 1
  integer(ik4) , parameter :: srf_xlat   = 2
  integer(ik4) , parameter :: srf_mask   = 3
  integer(ik4) , parameter :: srf_topo   = 4
  integer(ik4) , parameter :: srf_ps     = 5
  integer(ik4) , parameter :: srf_uvdrag = 6
  integer(ik4) , parameter :: srf_tg     = 7
  integer(ik4) , parameter :: srf_tlef   = 8
  integer(ik4) , parameter :: srf_tpr    = 9
  integer(ik4) , parameter :: srf_evp    = 10
  integer(ik4) , parameter :: srf_scv    = 11
  integer(ik4) , parameter :: srf_sena   = 12
  integer(ik4) , parameter :: srf_flw    = 13
  integer(ik4) , parameter :: srf_fsw    = 14
  integer(ik4) , parameter :: srf_fld    = 15
  integer(ik4) , parameter :: srf_sina   = 16
  integer(ik4) , parameter :: srf_prcv   = 17
  integer(ik4) , parameter :: srf_zpbl   = 18
  integer(ik4) , parameter :: srf_aldirs = 19
  integer(ik4) , parameter :: srf_aldifs = 20
  integer(ik4) , parameter :: srf_sund   = 21
  integer(ik4) , parameter :: srf_seaice = 22
  integer(ik4) , parameter :: srf_snowmelt = 23

  integer(ik4) , parameter :: srf_u10m   = 1
  integer(ik4) , parameter :: srf_v10m   = 2
  integer(ik4) , parameter :: srf_t2m    = 3
  integer(ik4) , parameter :: srf_q2m    = 4
  integer(ik4) , parameter :: srf_rh2m   = 5
  integer(ik4) , parameter :: srf_smw    = 6
  integer(ik4) , parameter :: srf_runoff = 7

  integer(ik4) , parameter :: sts_xlon   = 1
  integer(ik4) , parameter :: sts_xlat   = 2
  integer(ik4) , parameter :: sts_mask   = 3
  integer(ik4) , parameter :: sts_topo   = 4
  integer(ik4) , parameter :: sts_ps     = 5
  integer(ik4) , parameter :: sts_tgmax  = 6
  integer(ik4) , parameter :: sts_tgmin  = 7
  integer(ik4) , parameter :: sts_pcpmax = 8
  integer(ik4) , parameter :: sts_pcpavg = 9
  integer(ik4) , parameter :: sts_sund   = 10
  integer(ik4) , parameter :: sts_psmin  = 11
  integer(ik4) , parameter :: sts_psavg  = 12

  integer(ik4) , parameter :: sts_t2max  = 1
  integer(ik4) , parameter :: sts_t2min  = 2
  integer(ik4) , parameter :: sts_t2avg  = 3
  integer(ik4) , parameter :: sts_w10max = 4
  integer(ik4) , parameter :: sts_runoff = 5

  integer(ik4) , parameter :: sub_xlon   = 1
  integer(ik4) , parameter :: sub_xlat   = 2
  integer(ik4) , parameter :: sub_mask   = 3
  integer(ik4) , parameter :: sub_topo   = 4
  integer(ik4) , parameter :: sub_ps     = 5
  integer(ik4) , parameter :: sub_uvdrag = 6
  integer(ik4) , parameter :: sub_tg     = 7
  integer(ik4) , parameter :: sub_tlef   = 8
  integer(ik4) , parameter :: sub_evp    = 9
  integer(ik4) , parameter :: sub_scv    = 10
  integer(ik4) , parameter :: sub_sena   = 11

  integer(ik4) , parameter :: sub_u10m   = 1
  integer(ik4) , parameter :: sub_v10m   = 2
  integer(ik4) , parameter :: sub_t2m    = 3
  integer(ik4) , parameter :: sub_q2m    = 4
  integer(ik4) , parameter :: sub_smw    = 5
  integer(ik4) , parameter :: sub_runoff = 6

  integer(ik4) , parameter :: rad_xlon   = 1
  integer(ik4) , parameter :: rad_xlat   = 2
  integer(ik4) , parameter :: rad_mask   = 3
  integer(ik4) , parameter :: rad_topo   = 4
  integer(ik4) , parameter :: rad_ps     = 5
  integer(ik4) , parameter :: rad_frsa   = 6
  integer(ik4) , parameter :: rad_frla   = 7
  integer(ik4) , parameter :: rad_clrst  = 8
  integer(ik4) , parameter :: rad_clrss  = 9
  integer(ik4) , parameter :: rad_clrlt  = 10
  integer(ik4) , parameter :: rad_clrls  = 11
  integer(ik4) , parameter :: rad_solin  = 12
  integer(ik4) , parameter :: rad_sabtp  = 13
  integer(ik4) , parameter :: rad_totcf  = 14
  integer(ik4) , parameter :: rad_totcl  = 15
  integer(ik4) , parameter :: rad_totci  = 16
  integer(ik4) , parameter :: rad_firtp  = 17

  integer(ik4) , parameter :: rad_cld    = 1
  integer(ik4) , parameter :: rad_clwp   = 2
  integer(ik4) , parameter :: rad_qrs    = 3
  integer(ik4) , parameter :: rad_qrl    = 4

  integer(ik4) , parameter :: rad_taucl  = 1
  integer(ik4) , parameter :: rad_tauci  = 2

  integer(ik4) , parameter :: lak_xlon   = 1
  integer(ik4) , parameter :: lak_xlat   = 2
  integer(ik4) , parameter :: lak_mask   = 3
  integer(ik4) , parameter :: lak_topo   = 4
  integer(ik4) , parameter :: lak_ps     = 5
  integer(ik4) , parameter :: lak_tg     = 6
  integer(ik4) , parameter :: lak_tpr    = 7
  integer(ik4) , parameter :: lak_scv    = 8
  integer(ik4) , parameter :: lak_sena   = 9
  integer(ik4) , parameter :: lak_flw    = 10
  integer(ik4) , parameter :: lak_fsw    = 11
  integer(ik4) , parameter :: lak_fld    = 12
  integer(ik4) , parameter :: lak_sina   = 13
  integer(ik4) , parameter :: lak_aldirs = 14
  integer(ik4) , parameter :: lak_aldifs = 15
  integer(ik4) , parameter :: lak_evp    = 16
  integer(ik4) , parameter :: lak_ice    = 17

  integer(ik4) , parameter :: lak_tlake  = 1

  integer(ik4) , parameter :: opt_xlon     = 1
  integer(ik4) , parameter :: opt_xlat     = 2
  integer(ik4) , parameter :: opt_mask     = 3
  integer(ik4) , parameter :: opt_topo     = 4
  integer(ik4) , parameter :: opt_ps       = 5
  integer(ik4) , parameter :: opt_acstoarf = 6
  integer(ik4) , parameter :: opt_acstsrrf = 7
  integer(ik4) , parameter :: opt_acstalrf = 8
  integer(ik4) , parameter :: opt_acssrlrf = 9
  integer(ik4) , parameter :: opt_aod      = 10

  integer(ik4) , parameter :: opt_aext8    = 1
  integer(ik4) , parameter :: opt_assa8    = 2
  integer(ik4) , parameter :: opt_agfu8    = 3

  integer(ik4) , parameter :: che_xlon     = 1
  integer(ik4) , parameter :: che_xlat     = 2
  integer(ik4) , parameter :: che_mask     = 3
  integer(ik4) , parameter :: che_topo     = 4
  integer(ik4) , parameter :: che_ps       = 5
  integer(ik4) , parameter :: che_wdrflx   = 6
  integer(ik4) , parameter :: che_wdcflx   = 7
  integer(ik4) , parameter :: che_ddflx    = 8
  integer(ik4) , parameter :: che_emflx    = 9
  integer(ik4) , parameter :: che_ddvel    = 10
  integer(ik4) , parameter :: che_burden   = 11
  integer(ik4) , parameter :: che_pblten   = 12

  integer(ik4) , parameter :: che_mixrat   = 1
  integer(ik4) , parameter :: che_cheten   = 2
  integer(ik4) , parameter :: che_advhten  = 3
  integer(ik4) , parameter :: che_advvten  = 4
  integer(ik4) , parameter :: che_difhten  = 5
  integer(ik4) , parameter :: che_cuten    = 6
  integer(ik4) , parameter :: che_tuten    = 7
  integer(ik4) , parameter :: che_raiten   = 8
  integer(ik4) , parameter :: che_wasten   = 9
  integer(ik4) , parameter :: che_bdyten   = 10
  integer(ik4) , parameter :: che_sedten   = 11
  integer(ik4) , parameter :: che_emten    = 12

  integer(ik4) , parameter :: slab_xlon    = 1
  integer(ik4) , parameter :: slab_xlat    = 2
  integer(ik4) , parameter :: slab_mask    = 3
  integer(ik4) , parameter :: slab_topo    = 4
  integer(ik4) , parameter :: slab_ps      = 5

  integer(ik4) , public , parameter :: slab_qflx    = 1

  real(rk8) , pointer , dimension(:,:) :: io2d , io2dsg
  real(rk8) , pointer , dimension(:,:,:) :: io3d , io3dsg
  real(rk8) , pointer , dimension(:,:,:,:) :: io4d

  interface setup_var
    module procedure setup_var_2d
    module procedure setup_var_3d
    module procedure setup_var_4d
  end interface setup_var

  interface writevar_output_stream
    module procedure writevar2d_output_stream
    module procedure writevar3d_output_stream
  end interface

  contains

  subroutine init_output_streams(lparallel)
    implicit none
    logical , intent(in) :: lparallel
    integer(ik4) :: nstream , i , itr , vcount , kkz , n4dd
    type(varspan) :: vsize
    logical , dimension(natm2dvars) :: enable_atm2d_vars
    logical , dimension(natm3dvars) :: enable_atm3d_vars
    logical , dimension(nsrf2dvars) :: enable_srf2d_vars
    logical , dimension(nsrf3dvars) :: enable_srf3d_vars
    logical , dimension(nsts2dvars) :: enable_sts2d_vars
    logical , dimension(nsts3dvars) :: enable_sts3d_vars
    logical , dimension(nsub2dvars) :: enable_sub2d_vars
    logical , dimension(nsub3dvars) :: enable_sub3d_vars
    logical , dimension(nlak2dvars) :: enable_lak2d_vars
    logical , dimension(nlak3dvars) :: enable_lak3d_vars
    logical , dimension(nrad2dvars) :: enable_rad2d_vars
    logical , dimension(nrad3dvars) :: enable_rad3d_vars
    logical , dimension(nrad4dvars) :: enable_rad4d_vars
    logical , dimension(nopt2dvars) :: enable_opt2d_vars
    logical , dimension(nopt3dvars) :: enable_opt3d_vars
    logical , dimension(nche2dvars) :: enable_che2d_vars
    logical , dimension(nche3dvars) :: enable_che3d_vars

    parallel_out = lparallel

    nstream = 0

    if ( ifatm ) then
      nstream = nstream+1
      atm_stream = nstream
    end if
    if ( ifsrf ) then
      nstream = nstream+1
      srf_stream = nstream
    end if
    if ( ifsts ) then
      nstream = nstream+1
      sts_stream = nstream
    end if
    if ( ifsub ) then
      nstream = nstream+1
      sub_stream = nstream
    end if
    if ( iflak ) then
      nstream = nstream+1
      lak_stream = nstream
    end if
    if ( ifrad ) then
      nstream = nstream+1
      rad_stream = nstream
    end if
    if ( ifchem ) then
      nstream = nstream+1
      che_stream = nstream
      if ( iaerosol == 1 .and. ifopt ) then
        nstream = nstream+1
        opt_stream = nstream
      end if
    end if
    if ( ifslaboc ) then
      nstream = nstream+1
      slaboc_stream = nstream
    end if

    maxstreams = nstream
    allocate(outstream(maxstreams))

#ifdef DEBUG
    if ( myid == italk) &
      write(ndebug+myid,*) 'Enabling ',nstream,' output file streams'
#endif

    enabled_stream_loop: &
    do nstream = 1 , maxstreams

      vsize%j1 = jci1
      vsize%j2 = jci2
      vsize%i1 = ici1
      vsize%i2 = ici2
      vsize%k1 = 1
      vsize%n1 = 1
      vsize%n2 = 4

      if ( nstream == atm_stream ) then

        allocate(v2dvar_atm(natm2dvars))
        allocate(v3dvar_atm(natm3dvars))
        enable_atm2d_vars = enable_atm_vars(1:natm2dvars)
        enable_atm3d_vars = enable_atm_vars(natm2dvars+1:natmvars)

        ! This variables are always present

        call setup_common_vars(vsize,v2dvar_atm,atm_xlon, &
                  atm_xlat,atm_topo,atm_mask,atm_ps)

        ! The following may be enabled/disabled

        if ( enable_atm2d_vars(atm_tpr) ) then
          call setup_var(v2dvar_atm,atm_tpr,vsize,'pr','kg m-2 s-1', &
            'Total rain precipitation flux','precipitation_flux',.true., &
            'time: mean')
          atm_tpr_out => v2dvar_atm(atm_tpr)%rval
        end if
        if ( ipptls == 2 ) then
          if ( enable_atm2d_vars(atm_tsn) ) then
            call setup_var(v2dvar_atm,atm_tsn,vsize,'snw','kg m-2 s-1', &
              'Total snow precipitation flux','snow_flux',.true., &
              'time: mean')
            atm_tsn_out => v2dvar_atm(atm_tsn)%rval
          end if
        else
          enable_atm2d_vars(atm_tsn) = .false.
        end if
        if ( enable_atm2d_vars(atm_tgb) ) then
          call setup_var(v2dvar_atm,atm_tgb,vsize,'ts','K', &
            'Groud surface temperature','soil_temperature',.true.,'time: mean')
          atm_tgb_out => v2dvar_atm(atm_tgb)%rval
        end if
        if ( enable_atm2d_vars(atm_tsw) ) then
          call setup_var(v2dvar_atm,atm_tsw,vsize,'mrso','kg m-2', &
            'Total soil water','soil_moisture_content',.true., &
            'time: mean',l_fill=.true.)
          atm_tsw_out => v2dvar_atm(atm_tsw)%rval
        end if

        vsize%k2 = kz
        if ( enable_atm3d_vars(atm_u) ) then
          call setup_var(v3dvar_atm,atm_u,vsize,'ua','m s-1', &
            'Zonal component of wind (westerly)','eastward_wind',.true.)
          atm_u_out => v3dvar_atm(atm_u)%rval
        end if
        if ( enable_atm3d_vars(atm_v) ) then
          call setup_var(v3dvar_atm,atm_v,vsize,'va','m s-1', &
            'Meridional component of wind (southerly)','northward_wind',.true.)
          atm_v_out => v3dvar_atm(atm_v)%rval
        end if
        if ( enable_atm3d_vars(atm_t) ) then
          call setup_var(v3dvar_atm,atm_t,vsize,'ta','K', &
            'Air Temperature','air_temperature',.true.)
          atm_t_out => v3dvar_atm(atm_t)%rval
        end if
        if ( enable_atm3d_vars(atm_omega) ) then
          call setup_var(v3dvar_atm,atm_omega,vsize,'omega','hPa s-1', &
            'Pressure velocity','lagrangian_tendency_of_air_pressure',.true.)
          atm_omega_out => v3dvar_atm(atm_omega)%rval
        end if
        if ( enable_atm3d_vars(atm_qv) ) then
          call setup_var(v3dvar_atm,atm_qv,vsize,'qas','kg kg-1', &
            'Specific humidity in air','specific_humidity',.true.)
          atm_qv_out => v3dvar_atm(atm_qv)%rval
        end if
        if ( enable_atm3d_vars(atm_qc) ) then
          call setup_var(v3dvar_atm,atm_qc,vsize,'clw','kg kg-1', &
            'Mass fraction of cloud liquid water', &
            'mass_fraction_of_cloud_liquid_water_in_air',.true.)
          atm_qc_out => v3dvar_atm(atm_qc)%rval
        end if
        if ( enable_atm3d_vars(atm_rh) ) then
          call setup_var(v3dvar_atm,atm_rh,vsize,'rh','%', &
            'Relative Humidity', &
            'relative_humidity',.true.)
          atm_rh_out => v3dvar_atm(atm_rh)%rval
        end if
        if ( ipptls == 2 ) then
          if ( enable_atm3d_vars(atm_q_detr) ) then
            call setup_var(v3dvar_atm,atm_q_detr,vsize,'qdetr','kg/m^2', &
              'Detrainment', 'detrainment',.true.)
            atm_q_detr_out => v3dvar_atm(atm_q_detr)%rval
          end if
          if ( enable_atm3d_vars(atm_qr) ) then
            call setup_var(v3dvar_atm,atm_qr,vsize,'clr','kg kg-1', &
              'Mass fraction of rain', &
              'mass_fraction_of_rain',.true.)
            atm_qr_out => v3dvar_atm(atm_qr)%rval
          end if
          if ( enable_atm3d_vars(atm_qi) ) then
            call setup_var(v3dvar_atm,atm_qi,vsize,'cli','kg kg-1', &
              'Mass fraction of ice', &
              'mass_fraction_of_ice_in_air',.true.)
            atm_qi_out => v3dvar_atm(atm_qi)%rval
          end if
          if ( enable_atm3d_vars(atm_qs) ) then
            call setup_var(v3dvar_atm,atm_qs,vsize,'cls','kg kg-1', &
              'Mass fraction of snow', &
              'mass_fraction_of_snow_in_air',.true.)
            atm_qs_out => v3dvar_atm(atm_qs)%rval
          end if
          if ( enable_atm3d_vars(atm_zf) ) then
            call setup_var(v3dvar_atm,atm_zf,vsize,'zf','m', &
              'Height at full levels', &
              'height_full_levels',.true.)
            atm_zf_out => v3dvar_atm(atm_zf)%rval
          end if
          if ( enable_atm3d_vars(atm_zh) ) then
            call setup_var(v3dvar_atm,atm_zh,vsize,'zh','m', &
              'Height at half levels', &
              'height_half_levels',.true.)
            atm_zh_out => v3dvar_atm(atm_zh)%rval
          end if
          if ( enable_atm3d_vars(atm_pf) ) then
            call setup_var(v3dvar_atm,atm_pf,vsize,'pf','Pa', &
              'Pressure at full levels', &
              'pressure_full_levels',.true.)
            atm_pf_out => v3dvar_atm(atm_pf)%rval
          end if
          if ( enable_atm3d_vars(atm_ph) ) then
            call setup_var(v3dvar_atm,atm_ph,vsize,'ph','Pa', &
              'Pressure at half levels', &
              'pressure_half_levels',.true.)
            atm_ph_out => v3dvar_atm(atm_ph)%rval
          end if
          if ( idiag > 0) then
            if ( enable_atm3d_vars(atm_rainls) ) then
              call setup_var(v3dvar_atm,atm_rainls,vsize, &
                 'rainls','kg m-2 s-1', &
                 'Large scale precipitation at each level', &
                 'large_scale_precipitation',.true.)
              atm_rainls_out => v3dvar_atm(atm_rainls)%rval
            end if
            if ( enable_atm3d_vars(atm_raincc) ) then
              call setup_var(v3dvar_atm,atm_raincc,vsize, &
                'raincc','kg m-2 s-1', &
                'Convective precipitation at each level', &
                'convective_precipitation',.true.)
              atm_raincc_out => v3dvar_atm(atm_raincc)%rval
            end if
          else
            enable_atm3d_vars(atm_rainls) = .false.
            enable_atm3d_vars(atm_raincc) = .false.
          end if
          if (stats) then
            ! stats variables
            if ( enable_atm3d_vars(atm_stats_supw) ) then
              call setup_var(v3dvar_atm,atm_stats_supw,vsize,'st_supw','', &
              '',&
              '',.true.)
              atm_stats_supw_out => v3dvar_atm(atm_stats_supw)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_supc) ) then
              call setup_var(v3dvar_atm,atm_stats_supc,vsize,'st_supc','', &
              '',&
              '',.true.)
              atm_stats_supc_out => v3dvar_atm(atm_stats_supc)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_detw) ) then
              call setup_var(v3dvar_atm,atm_stats_detw,vsize,'st_detw','', &
              '',&
              '',.true.)
              atm_stats_detw_out => v3dvar_atm(atm_stats_detw)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_detc) ) then
              call setup_var(v3dvar_atm,atm_stats_detc,vsize,'st_detc','', &
              '',&
              '',.true.)
              atm_stats_detc_out => v3dvar_atm(atm_stats_detc)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_erow) ) then
              call setup_var(v3dvar_atm,atm_stats_erow,vsize,'st_erow','', &
              '',&
              '',.true.)
              atm_stats_erow_out => v3dvar_atm(atm_stats_erow)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_eroc) ) then
              call setup_var(v3dvar_atm,atm_stats_eroc,vsize,'st_eroc','', &
              '',&
              '',.true.)
              atm_stats_eroc_out => v3dvar_atm(atm_stats_eroc)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_evw) )  then
              call setup_var(v3dvar_atm,atm_stats_evw,vsize,'st_evw','', &
              '',&
              '',.true.)
              atm_stats_evw_out => v3dvar_atm(atm_stats_evw)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_evc) ) then
              call setup_var(v3dvar_atm,atm_stats_evc,vsize,'st_evc','', &
              '',&
              '',.true.)
              atm_stats_evc_out => v3dvar_atm(atm_stats_evc)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_con1w) ) then
              call setup_var(v3dvar_atm,atm_stats_con1w,vsize,'st_con1w','', &
              '',&
              '',.true.)
              atm_stats_con1w_out => v3dvar_atm(atm_stats_con1w)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_con1c) ) then
              call setup_var(v3dvar_atm,atm_stats_con1c,vsize,'st_con1c','', &
              '',&
              '',.true.)
              atm_stats_con1c_out => v3dvar_atm(atm_stats_con1c)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_con2w) ) then
              call setup_var(v3dvar_atm,atm_stats_con2w,vsize,'st_con2w','', &
              '',&
              '',.true.)
             atm_stats_con2w_out => v3dvar_atm(atm_stats_con2w)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_con2c) ) then
              call setup_var(v3dvar_atm,atm_stats_con2c,vsize,'st_con2c','', &
              '',&
              '',.true.)
             atm_stats_con2c_out => v3dvar_atm(atm_stats_con2c)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_dep) ) then
              call setup_var(v3dvar_atm,atm_stats_dep,vsize,'st_dep','', &
              '',&
              '',.true.)
              atm_stats_dep_out => v3dvar_atm(atm_stats_dep)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_melt) ) then
              call setup_var(v3dvar_atm,atm_stats_melt,vsize,'st_mlt','', &
              '',&
              '',.true.)
              atm_stats_melt_out=> v3dvar_atm(atm_stats_melt)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_frz) ) then
              call setup_var(v3dvar_atm,atm_stats_frz,vsize,'st_frz','', &
              '',&
              '',.true.)
              atm_stats_frz_out=> v3dvar_atm(atm_stats_frz)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_rainev) ) then
              call setup_var(v3dvar_atm,atm_stats_rainev,vsize,'st_ev_rn','', &
              '',&
              '',.true.)
              atm_stats_rainev_out=> v3dvar_atm(atm_stats_rainev)%rval
            end if
            if ( enable_atm3d_vars(atm_stats_snowev) ) then
              call setup_var(v3dvar_atm,atm_stats_snowev,vsize,'st_ev_sn','', &
              '',&
              '',.true.)
              atm_stats_snowev_out=> v3dvar_atm(atm_stats_snowev)%rval
            end if
          else!stats
            enable_atm3d_vars(atm_stats_supw:atm_stats_snowev) = .false.
          end if
        else
          enable_atm3d_vars(atm_q_detr) = .false.
          enable_atm3d_vars(atm_qr:atm_qs) = .false.
          enable_atm3d_vars(atm_zf) = .false.
          enable_atm3d_vars(atm_zh) = .false.
          enable_atm3d_vars(atm_pf) = .false.
          enable_atm3d_vars(atm_ph) = .false.
          enable_atm3d_vars(atm_rainls) = .false.
          enable_atm3d_vars(atm_raincc) = .false.
          enable_atm3d_vars(atm_stats_supw:atm_stats_snowev) = .false.
        end if

        if ( ibltyp == 2 .or. ibltyp == 99 ) then
          if ( enable_atm3d_vars(atm_tke) ) then
            call setup_var(v3dvar_atm,atm_tke,vsize,'tke','m2 s2', &
              'Turbulent Kinetic Energy','turbulent_kinetic_energy', .true.)
            atm_tke_out => v3dvar_atm(atm_tke)%rval
          end if
          if ( enable_atm3d_vars(atm_kth) ) then
            call setup_var(v3dvar_atm,atm_kth,vsize,'kth','m2 s-1', &
              'Vertical Turbulent Viscosity','vertical_momentum_diffusivity', &
              .true.)
            atm_kth_out => v3dvar_atm(atm_kth)%rval
          end if
          if ( enable_atm3d_vars(atm_kzm) ) then
            call setup_var(v3dvar_atm,atm_kzm,vsize,'kzm','m2 s-1', &
              'Vertical Turbulent Diffusivity','vertical_scalar_diffusivity', &
              .true.)
            atm_kzm_out => v3dvar_atm(atm_kzm)%rval
          end if
        else
          enable_atm3d_vars(atm_tke:atm_kzm) = .false.
        end if
        if ( idiag > 0 ) then
          ! FAB : flag properly
          if ( enable_atm3d_vars(atm_tten_adh) ) then
            call setup_var(v3dvar_atm,atm_tten_adh,vsize,'ttenadh','K.s-1', &
             'temperature_tendency_due_to_horizontal_advection', &
             'Temperature tendency due to horizontal advection',.true.)
            atm_tten_adh_out => v3dvar_atm(atm_tten_adh)%rval
          end if
          if ( enable_atm3d_vars(atm_tten_adv) ) then
            call setup_var(v3dvar_atm,atm_tten_adv,vsize,'ttenadv','K.s-1', &
             'temperature_tendency_due_to_vertical_advection', &
             'Temperature tendency due to vertical advection',.true.)
            atm_tten_adv_out => v3dvar_atm(atm_tten_adv)%rval
          end if
          if ( enable_atm3d_vars(atm_tten_tbl) ) then
            call setup_var(v3dvar_atm,atm_tten_tbl,vsize,'ttentbl','K.s-1', &
             'temperature_tendency_due_to_surface_boundary_layer', &
             'Temperature tendency due to surface boundary layer',.true.)
            atm_tten_tbl_out => v3dvar_atm(atm_tten_tbl)%rval
          end if
          if ( enable_atm3d_vars(atm_tten_dif) ) then
            call setup_var(v3dvar_atm,atm_tten_dif,vsize,'ttendif','K.s-1', &
             'temperature_tendency_due_to_diffusion', &
             'Temperature tendency due to diffusion',.true.)
            atm_tten_dif_out => v3dvar_atm(atm_tten_dif)%rval
          end if
          if ( enable_atm3d_vars(atm_tten_bdy) ) then
            call setup_var(v3dvar_atm,atm_tten_bdy,vsize,'ttenbdy','K.s-1', &
             'temperature_tendency_due_to_boundary_conditions', &
             'Temperature tendency due to boundary conditions',.true.)
            atm_tten_bdy_out => v3dvar_atm(atm_tten_bdy)%rval
          end if
          if ( enable_atm3d_vars(atm_tten_con) ) then
            call setup_var(v3dvar_atm,atm_tten_con,vsize,'ttencon','K.s-1', &
             'temperature_tendency_due_to_convection', &
             'Temperature tendency due to convection',.true.)
            atm_tten_con_out => v3dvar_atm(atm_tten_con)%rval
          end if
          if ( enable_atm3d_vars(atm_tten_adi) ) then
            call setup_var(v3dvar_atm,atm_tten_adi,vsize,'ttenadi','K.s-1', &
             'temperature_tendency_due_to_adiabatic', &
             'Temperature tendency due to adiabatic',.true.)
            atm_tten_adi_out => v3dvar_atm(atm_tten_adi)%rval
          end if
          if ( enable_atm3d_vars(atm_tten_rad) ) then
            call setup_var(v3dvar_atm,atm_tten_rad,vsize,'ttenrad','K.s-1', &
             'temperature_tendency_due_to_radiation_heating', &
             'Temperature tendency due to radiation heating',.true.)
            atm_tten_rad_out => v3dvar_atm(atm_tten_rad)%rval
          end if
          if ( enable_atm3d_vars(atm_tten_lsc) ) then
            call setup_var(v3dvar_atm,atm_tten_lsc,vsize,'ttenlsc','K.s-1', &
             'temperature_tendency_due_to_large_scale_latent_heat_exchange', &
             'Temperature tendency due to large scale latent heat exchange', &
             .true.)
            atm_tten_lsc_out => v3dvar_atm(atm_tten_lsc)%rval
          end if
          if ( enable_atm3d_vars(atm_qten_adh) ) then
            call setup_var(v3dvar_atm,atm_qten_adh,vsize,'qtenadh','s-1', &
             'mixing_ratio_tendency_due_to_horizontal_advection', &
             'Mixing ratio tendency due to horizontal advection',.true.)
            atm_qten_adh_out => v3dvar_atm(atm_qten_adh)%rval
          end if
          if ( enable_atm3d_vars(atm_qten_adv) ) then
            call setup_var(v3dvar_atm,atm_qten_adv,vsize,'qtenadv','s-1', &
             'mixing_ratio_tendency_due_to_vertical_advection', &
             'Mixing ratio tendency due to vertical advection',.true.)
            atm_qten_adv_out => v3dvar_atm(atm_qten_adv)%rval
          end if
          if ( enable_atm3d_vars(atm_qten_tbl) ) then
            call setup_var(v3dvar_atm,atm_qten_tbl,vsize,'qtentbl','s-1', &
             'mixing_ratio_tendency_due_to_surface_boundary_layer', &
             'Mixing ratio tendency due to surface boundary layer',.true.)
            atm_qten_tbl_out => v3dvar_atm(atm_qten_tbl)%rval
          end if
          if ( enable_atm3d_vars(atm_qten_dif) ) then
            call setup_var(v3dvar_atm,atm_qten_dif,vsize,'qtendif','s-1', &
             'mixing_ratio_tendency_due_to_diffusion', &
             'Mixing ratio tendency due to diffusion',.true.)
            atm_qten_dif_out => v3dvar_atm(atm_qten_dif)%rval
          end if
          if ( enable_atm3d_vars(atm_qten_bdy) ) then
            call setup_var(v3dvar_atm,atm_qten_bdy,vsize,'qtenbdy','s-1', &
             'mixing_ratio_tendency_due_to_boundary_conditions', &
             'Mixing ratio tendency due to boundary conditions',.true.)
            atm_qten_bdy_out => v3dvar_atm(atm_qten_bdy)%rval
          end if
          if ( enable_atm3d_vars(atm_qten_con) ) then
            call setup_var(v3dvar_atm,atm_qten_con,vsize,'qtencon','s-1', &
             'mixing_ratio_tendency_due_to_convection', &
             'Mixing ratio tendency due to convection',.true.)
            atm_qten_con_out => v3dvar_atm(atm_qten_con)%rval
          end if
          if ( enable_atm3d_vars(atm_qten_adi) ) then
            call setup_var(v3dvar_atm,atm_qten_adi,vsize,'qtenadi','s-1', &
             'mixing_ratio_tendency_due_to_adiabatic', &
             'Mixing ratio tendency due to adiabatic',.true.)
            atm_qten_adi_out => v3dvar_atm(atm_qten_adi)%rval
          end if
          if ( enable_atm3d_vars(atm_qten_rad) ) then
            call setup_var(v3dvar_atm,atm_qten_rad,vsize,'qtenrad','s-1', &
             'mixing_ratio_tendency_due_to_radiation_heating', &
             'Mixing ratio tendency due to radiation heating',.true.)
            atm_qten_rad_out => v3dvar_atm(atm_qten_rad)%rval
          end if
          if ( enable_atm3d_vars(atm_qten_lsc) ) then
            call setup_var(v3dvar_atm,atm_qten_lsc,vsize,'qtenlsc','s-1', &
             'mixing_ratio_tendency_due_to_large_scale_latent_heat_exchange', &
             'Mixing ratio tendency due to large scale latent heat exchange', &
             .true.)
            atm_qten_lsc_out => v3dvar_atm(atm_qten_lsc)%rval
          end if
        else
          enable_atm3d_vars(atm_tten_adh:atm_qten_lsc) = .false.
        end if

        enable_atm_vars(1:natm2dvars) = enable_atm2d_vars
        enable_atm_vars(natm2dvars+1:natmvars) = enable_atm3d_vars
        outstream(atm_stream)%nvar = countvars(enable_atm_vars,natmvars)
        allocate(outstream(atm_stream)%ncvars%vlist(outstream(atm_stream)%nvar))
        outstream(atm_stream)%nfiles = 1
        allocate(outstream(atm_stream)%ncout(outstream(atm_stream)%nfiles))
        allocate(outstream(atm_stream)%cname_base(outstream(atm_stream)%nfiles))
        outstream(atm_stream)%cname_base(1) = 'ATM'

        vcount = 1
        do i = 1 , natm2dvars
          if ( enable_atm_vars(i) ) then
            outstream(atm_stream)%ncvars%vlist(vcount)%vp => v2dvar_atm(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , natm3dvars
          if ( enable_atm_vars(i+natm2dvars) ) then
            outstream(atm_stream)%ncvars%vlist(vcount)%vp => v3dvar_atm(i)
            vcount = vcount + 1
          end if
        end do
        outstream(atm_stream)%jl1 = vsize%j1
        outstream(atm_stream)%jl2 = vsize%j2
        outstream(atm_stream)%il1 = vsize%i1
        outstream(atm_stream)%il2 = vsize%i2
        outstream(atm_stream)%jg1 = jout1
        outstream(atm_stream)%jg2 = jout2
        outstream(atm_stream)%ig1 = iout1
        outstream(atm_stream)%ig2 = iout2
      end if

      if ( nstream == srf_stream ) then

        allocate(v2dvar_srf(nsrf2dvars))
        allocate(v3dvar_srf(nsrf3dvars))
        enable_srf2d_vars = enable_srf_vars(1:nsrf2dvars)
        enable_srf3d_vars = enable_srf_vars(nsrf2dvars+1:nsrfvars)

        ! This variables are always present

        call setup_common_vars(vsize,v2dvar_srf,srf_xlon, &
                  srf_xlat,srf_topo,srf_mask,srf_ps)

        ! The following may be enabled/disabled

        if ( enable_srf2d_vars(srf_uvdrag) ) then
          call setup_var(v2dvar_srf,srf_uvdrag,vsize,'tau','N m-2', &
            'Surface wind stress','surface_downward_stress',.true.)
          srf_uvdrag_out => v2dvar_srf(srf_uvdrag)%rval
        end if
        if ( enable_srf2d_vars(srf_tg) ) then
          call setup_var(v2dvar_srf,srf_tg,vsize,'ts','K', &
            'Ground surface temperature','surface_temperature',.true.)
          srf_tg_out => v2dvar_srf(srf_tg)%rval
        end if
        if ( enable_srf2d_vars(srf_tlef) ) then
          call setup_var(v2dvar_srf,srf_tlef,vsize,'tf','K', &
            'Foliage canopy temperature','canopy_temperature',  &
            .true.,l_fill=.true.)
          srf_tlef_out => v2dvar_srf(srf_tlef)%rval
        end if
        if ( enable_srf2d_vars(srf_tpr) ) then
          call setup_var(v2dvar_srf,srf_tpr,vsize,'pr','kg m-2 s-1', &
            'Total precipitation flux','precipitation_flux',.true., &
            'time: mean')
          srf_tpr_out => v2dvar_srf(srf_tpr)%rval
        end if
        if ( enable_srf2d_vars(srf_evp) ) then
          call setup_var(v2dvar_srf,srf_evp,vsize,'evspsbl','kg m-2 s-1', &
            'Total evapotranspiration flux','water_evaporation_flux',.true., &
            'time: mean')
          srf_evp_out => v2dvar_srf(srf_evp)%rval
        end if
        if ( enable_srf2d_vars(srf_scv) ) then
          call setup_var(v2dvar_srf,srf_scv,vsize,'snv','kg m-2', &
            'Liquid water equivalent of snow thickness', &
            'lwe_thickness_of_surface_snow_amount',.true.,'time: mean', &
            l_fill=.true.)
          srf_scv_out => v2dvar_srf(srf_scv)%rval
        end if
        if ( enable_srf2d_vars(srf_sena) ) then
          call setup_var(v2dvar_srf,srf_sena,vsize,'hfss','W m-2', &
            'Sensible heat flux','surface_upward_sensible_heat_flux', &
            .true.,'time: mean')
          srf_sena_out => v2dvar_srf(srf_sena)%rval
        end if
        if ( enable_srf2d_vars(srf_flw) ) then
          call setup_var(v2dvar_srf,srf_flw,vsize,'rsnl','W m-2', &
            'Net upward longwave energy flux', &
            'net_upward_longwave_flux_in_air',.true.,'time: mean')
          srf_flw_out => v2dvar_srf(srf_flw)%rval
        end if
        if ( enable_srf2d_vars(srf_fsw) ) then
          call setup_var(v2dvar_srf,srf_fsw,vsize,'rsns','W m-2', &
            'Net downward shortwave energy flux', &
            'net_downward_shortwave_flux_in_air',.true.,'time: mean')
          srf_fsw_out => v2dvar_srf(srf_fsw)%rval
        end if
        if ( enable_srf2d_vars(srf_fld) ) then
          call setup_var(v2dvar_srf,srf_fld,vsize,'rsdl','W m-2', &
            'Surface downward longwave flux in air', &
            'surface_downwelling_longwave_flux_in_air',.true.,'time: mean')
          srf_fld_out => v2dvar_srf(srf_fld)%rval
        end if
        if ( enable_srf2d_vars(srf_sina) ) then
          call setup_var(v2dvar_srf,srf_sina,vsize,'rsds','W m-2', &
            'Surface downward shortwave flux in air', &
            'surface_downwelling_shortwave_flux_in_air',.true.,'time: mean')
          srf_sina_out => v2dvar_srf(srf_sina)%rval
        end if
        if ( enable_srf2d_vars(srf_prcv) ) then
          call setup_var(v2dvar_srf,srf_prcv,vsize,'prc','kg m-2 s-1', &
            'Convective precipitation flux','convective_rainfall_flux', &
            .true.,'time: mean')
          srf_prcv_out => v2dvar_srf(srf_prcv)%rval
        end if
        if ( enable_srf2d_vars(srf_zpbl) ) then
          call setup_var(v2dvar_srf,srf_zpbl,vsize,'zmla','m', &
            'Atmospheric Boundary Layer thickness', &
            'atmosphere_boundary_layer_thickness',.true.)
          srf_zpbl_out => v2dvar_srf(srf_zpbl)%rval
        end if
        if ( enable_srf2d_vars(srf_aldirs) ) then
          call setup_var(v2dvar_srf,srf_aldirs,vsize,'aldirs','1', &
            'Surface albedo to direct shortwave radiation', &
            'surface_albedo_short_wave_direct',.true.)
          srf_aldirs_out => v2dvar_srf(srf_aldirs)%rval
        end if
        if ( enable_srf2d_vars(srf_aldifs) ) then
          call setup_var(v2dvar_srf,srf_aldifs,vsize,'aldifs','1', &
            'Surface albedo to diffuse shortwave radiation', &
            'surface_albedo_short_wave_diffuse',.true.)
          srf_aldifs_out => v2dvar_srf(srf_aldifs)%rval
        end if
        if ( enable_srf2d_vars(srf_sund) ) then
          call setup_var(v2dvar_srf,srf_sund,vsize,'sund','s', &
            'Duration of sunshine','duration_of_sunshine',.true.,'time: sum')
          srf_sund_out => v2dvar_srf(srf_sund)%rval
        end if
        if ( iseaice == 1 ) then
          if ( enable_srf2d_vars(srf_seaice) ) then
            call setup_var(v2dvar_srf,srf_seaice,vsize,'seaice','m', &
              'Sea ice cover','seaice_depth',.true.)
            srf_seaice_out => v2dvar_srf(srf_seaice)%rval
          end if
        else
          enable_srf2d_vars(srf_seaice) = .false.
        end if
        if ( enable_srf2d_vars(srf_snowmelt) ) then
          call setup_var(v2dvar_srf,srf_snowmelt,vsize,'smelt','kg m-2', &
            'Snow Melt','surface_snow_melt_amount',.true.,'time: mean')
          srf_snowmelt_out => v2dvar_srf(srf_snowmelt)%rval
        end if

        vsize%k2 = 1
        v3dvar_srf(srf_u10m)%axis = 'xyw'
        v3dvar_srf(srf_v10m)%axis = 'xyw'
        v3dvar_srf(srf_t2m)%axis = 'xy2'
        v3dvar_srf(srf_q2m)%axis = 'xy2'
        v3dvar_srf(srf_rh2m)%axis = 'xy2'
        if ( enable_srf3d_vars(srf_u10m) ) then
          call setup_var(v3dvar_srf,srf_u10m,vsize,'uas','m s-1', &
            'Anemometric zonal (westerly) wind component', &
            'eastward_wind',.true.)
          srf_u10m_out => v3dvar_srf(srf_u10m)%rval
        end if
        if ( enable_srf3d_vars(srf_v10m) ) then
          call setup_var(v3dvar_srf,srf_v10m,vsize,'vas','m s-1', &
            'Anenometric meridional (southerly) wind component' , &
            'northward_wind',.true.)
          srf_v10m_out => v3dvar_srf(srf_v10m)%rval
        end if
        if ( enable_srf3d_vars(srf_t2m) ) then
          call setup_var(v3dvar_srf,srf_t2m,vsize,'tas','K', &
            'Near surface air temperature','air_temperature',.true.)
          srf_t2m_out => v3dvar_srf(srf_t2m)%rval
        end if
        if ( enable_srf3d_vars(srf_q2m) ) then
          call setup_var(v3dvar_srf,srf_q2m,vsize,'qas','1', &
            'Near surface air specific humidity','specific_humidity',.true.)
          srf_q2m_out => v3dvar_srf(srf_q2m)%rval
        end if
        if ( enable_srf3d_vars(srf_rh2m) ) then
          call setup_var(v3dvar_srf,srf_rh2m,vsize,'hurs','%', &
            'Near surface relative humidity','relative_humidity',.true.)
          srf_rh2m_out => v3dvar_srf(srf_rh2m)%rval
        end if

        vsize%k2 = 2
        v3dvar_srf(srf_smw)%axis = 'xys'
        v3dvar_srf(srf_runoff)%axis = 'xys'
        if ( enable_srf3d_vars(srf_smw) ) then
          call setup_var(v3dvar_srf,srf_smw,vsize,'mrso','kg m-2', &
            'Moisture content of the soil layers', &
            'soil_moisture_content_in_layers',.true.,l_fill=.true.)
          srf_smw_out => v3dvar_srf(srf_smw)%rval
        end if
        if ( enable_srf3d_vars(srf_runoff) ) then
          call setup_var(v3dvar_srf,srf_runoff,vsize,'mrro','kg m-2 s-1', &
            'Runoff flux','runoff_flux',.true.,'time: mean',l_fill=.true.)
          srf_runoff_out => v3dvar_srf(srf_runoff)%rval
        end if

        enable_srf_vars(1:nsrf2dvars) = enable_srf2d_vars
        enable_srf_vars(nsrf2dvars+1:nsrfvars) = enable_srf3d_vars
        outstream(srf_stream)%nvar = countvars(enable_srf_vars,nsrfvars)
        allocate(outstream(srf_stream)%ncvars%vlist(outstream(srf_stream)%nvar))
        outstream(srf_stream)%nfiles = 1
        allocate(outstream(srf_stream)%ncout(outstream(srf_stream)%nfiles))
        allocate(outstream(srf_stream)%cname_base(outstream(srf_stream)%nfiles))
        outstream(srf_stream)%cname_base(1) = 'SRF'

        vcount = 1
        do i = 1 , nsrf2dvars
          if ( enable_srf_vars(i) ) then
            outstream(srf_stream)%ncvars%vlist(vcount)%vp => v2dvar_srf(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nsrf3dvars
          if ( enable_srf_vars(i+nsrf2dvars) ) then
            outstream(srf_stream)%ncvars%vlist(vcount)%vp => v3dvar_srf(i)
            vcount = vcount + 1
          end if
        end do
        outstream(srf_stream)%jl1 = vsize%j1
        outstream(srf_stream)%jl2 = vsize%j2
        outstream(srf_stream)%il1 = vsize%i1
        outstream(srf_stream)%il2 = vsize%i2
        outstream(srf_stream)%jg1 = jout1
        outstream(srf_stream)%jg2 = jout2
        outstream(srf_stream)%ig1 = iout1
        outstream(srf_stream)%ig2 = iout2
      end if

      if ( nstream == sts_stream ) then

        allocate(v2dvar_sts(nsts2dvars))
        allocate(v3dvar_sts(nsts3dvars))
        enable_sts2d_vars = enable_sts_vars(1:nsts2dvars)
        enable_sts3d_vars = enable_sts_vars(nsts2dvars+1:nstsvars)

        ! This variables are always present

        call setup_common_vars(vsize,v2dvar_sts,sts_xlon, &
                  sts_xlat,sts_topo,sts_mask,sts_ps)

        ! The following may be enabled/disabled

        if ( enable_sts2d_vars(sts_tgmax) ) then
          call setup_var(v2dvar_sts,sts_tgmax,vsize,'tsmax','K', &
            'Maximum surface temperature','surface_temperature', &
            .true.,'time: maximum')
          sts_tgmax_out => v2dvar_sts(sts_tgmax)%rval
        end if
        if ( enable_sts2d_vars(sts_tgmin) ) then
          call setup_var(v2dvar_sts,sts_tgmin,vsize,'tsmin','K', &
            'Minimum surface temperature','surface_temperature', &
            .true.,'time: minimum')
          sts_tgmin_out => v2dvar_sts(sts_tgmin)%rval
        end if
        if ( enable_sts2d_vars(sts_pcpmax) ) then
          call setup_var(v2dvar_sts,sts_pcpmax,vsize,'prmax','kg m-2 s-1', &
            'Maximum total precipitation flux','precipitation_flux', &
            .true.,'time: maximum')
          sts_pcpmax_out => v2dvar_sts(sts_pcpmax)%rval
        end if
        if ( enable_sts2d_vars(sts_pcpavg) ) then
          call setup_var(v2dvar_sts,sts_pcpavg,vsize,'pr','kg m-2 s-1', &
            'Mean total precipitation flux','precipitation_flux', &
            .true.,'time: mean')
          sts_pcpavg_out => v2dvar_sts(sts_pcpavg)%rval
        end if
        if ( enable_sts2d_vars(sts_sund) ) then
          call setup_var(v2dvar_sts,sts_sund,vsize,'sund','s', &
            'Duration of sunshine','duration_of_sunshine',.true.,'time: sum')
          sts_sund_out => v2dvar_sts(sts_sund)%rval
        end if
        if ( enable_sts2d_vars(sts_psmin) ) then
          call setup_var(v2dvar_sts,sts_psmin,vsize,'psmin','hPa', &
            'Minimum of surface pressure','air_pressure',.true.,'time: minimum')
          sts_psmin_out => v2dvar_sts(sts_psmin)%rval
        end if
        if ( enable_sts2d_vars(sts_psavg) ) then
          call setup_var(v2dvar_sts,sts_psavg,vsize,'psavg','hPa', &
            'Mean surface pressure','air_pressure',.true.,'time: mean')
          sts_psavg_out => v2dvar_sts(sts_psavg)%rval
        end if

        vsize%k2 = 1
        v3dvar_sts(sts_t2max)%axis = 'xy2'
        v3dvar_sts(sts_t2min)%axis = 'xy2'
        v3dvar_sts(sts_t2avg)%axis = 'xy2'
        v3dvar_sts(sts_w10max)%axis = 'xyw'
        if ( enable_sts3d_vars(sts_t2max) ) then
          call setup_var(v3dvar_sts,sts_t2max,vsize,'tasmax','K', &
            'Maximum 2 meter temperature','air_temperature',.true., &
            'time: maximum')
          sts_t2max_out => v3dvar_sts(sts_t2max)%rval
        end if
        if ( enable_sts3d_vars(sts_t2min) ) then
          call setup_var(v3dvar_sts,sts_t2min,vsize,'tasmin','K', &
            'Minimum 2 meter temperature','air_temperature',.true., &
            'time: minimum')
          sts_t2min_out => v3dvar_sts(sts_t2min)%rval
        end if
        if ( enable_sts3d_vars(sts_t2avg) ) then
          call setup_var(v3dvar_sts,sts_t2avg,vsize,'tas','K', &
            'Mean 2 meter temperature','air_temperature',.true., &
            'time: mean')
          sts_t2avg_out => v3dvar_sts(sts_t2avg)%rval
        end if
        if ( enable_sts3d_vars(sts_w10max) ) then
          call setup_var(v3dvar_sts,sts_w10max,vsize,'sfcWindmax','m s-1', &
            'Maximum speed of 10m wind','wind_speed',.true.,'time: maximum')
          sts_w10max_out => v3dvar_sts(sts_w10max)%rval
        end if
        vsize%k2 = 2
        v3dvar_sts(sts_runoff)%axis = 'xys'
        if ( enable_sts3d_vars(sts_runoff) ) then
          call setup_var(v3dvar_sts,sts_runoff,vsize,'mrro','kg m-2 s-1', &
            'Runoff flux','runoff_flux',.true.,'time: mean',l_fill=.true.)
          sts_runoff_out => v3dvar_sts(sts_runoff)%rval
        end if

        enable_sts_vars(1:nsts2dvars) = enable_sts2d_vars
        enable_sts_vars(nsts2dvars+1:nstsvars) = enable_sts3d_vars
        outstream(sts_stream)%nvar = countvars(enable_sts_vars,nstsvars)
        allocate(outstream(sts_stream)%ncvars%vlist(outstream(sts_stream)%nvar))
        outstream(sts_stream)%nfiles = 1
        allocate(outstream(sts_stream)%ncout(outstream(sts_stream)%nfiles))
        allocate(outstream(sts_stream)%cname_base(outstream(sts_stream)%nfiles))
        outstream(sts_stream)%cname_base(1) = 'STS'

        vcount = 1
        do i = 1 , nsts2dvars
          if ( enable_sts_vars(i) ) then
            outstream(sts_stream)%ncvars%vlist(vcount)%vp => v2dvar_sts(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nsts3dvars
          if ( enable_sts_vars(i+nsts2dvars) ) then
            outstream(sts_stream)%ncvars%vlist(vcount)%vp => v3dvar_sts(i)
            vcount = vcount + 1
          end if
        end do
        outstream(sts_stream)%jl1 = vsize%j1
        outstream(sts_stream)%jl2 = vsize%j2
        outstream(sts_stream)%il1 = vsize%i1
        outstream(sts_stream)%il2 = vsize%i2
        outstream(sts_stream)%jg1 = jout1
        outstream(sts_stream)%jg2 = jout2
        outstream(sts_stream)%ig1 = iout1
        outstream(sts_stream)%ig2 = iout2
      end if

      if ( nstream == sub_stream ) then

        allocate(v2dvar_sub(nsub2dvars))
        allocate(v3dvar_sub(nsub3dvars))
        enable_sub2d_vars = enable_sub_vars(1:nsub2dvars)
        enable_sub3d_vars = enable_sub_vars(nsub2dvars+1:nsubvars)

        outstream(nstream)%opar%l_subgrid = .true.
        outstream(nstream)%l_sub = .true.

        vsize%j1 = (vsize%j1-1)*nsg+1
        vsize%j2 = vsize%j2*nsg
        vsize%i1 = (vsize%i1-1)*nsg+1
        vsize%i2 = vsize%i2*nsg

        ! This variables are always present

        call setup_var(v2dvar_sub,sub_xlon,vsize,'xlon','degrees_east', &
          'Longitude on Cross Points','longitude')
        call setup_var(v2dvar_sub,sub_xlat,vsize,'xlat','degrees_north', &
          'Latitude on Cross Points','latitude')
        call setup_var(v2dvar_sub,sub_mask,vsize,'mask','1', &
          'Land Mask','land_binary_mask')
        call setup_var(v2dvar_sub,sub_topo,vsize,'topo','m', &
          'Surface Model Elevation','surface_altitude')
        call setup_var(v2dvar_sub,sub_ps,vsize,'ps','hPa', &
          'Surface Pressure','surface_air_pressure',.true.,'time: mean')

        sub_xlon_out => v2dvar_sub(sub_xlon)%rval
        sub_xlat_out => v2dvar_sub(sub_xlat)%rval
        sub_mask_out => v2dvar_sub(sub_mask)%rval
        sub_topo_out => v2dvar_sub(sub_topo)%rval
        sub_ps_out => v2dvar_sub(sub_ps)%rval

        ! The following may be enabled/disabled

        if ( enable_sub2d_vars(sub_uvdrag) ) then
          call setup_var(v2dvar_sub,sub_uvdrag,vsize,'tau','N m-2', &
            'Surface wind stress','surface_downward_stress',.true.)
          sub_uvdrag_out => v2dvar_sub(sub_uvdrag)%rval
        end if
        if ( enable_sub2d_vars(sub_tg) ) then
          call setup_var(v2dvar_sub,sub_tg,vsize,'ts','K', &
            'Ground temperature','surface_temperature',.true.)
          sub_tg_out => v2dvar_sub(sub_tg)%rval
        end if
        if ( enable_sub2d_vars(sub_tlef) ) then
          call setup_var(v2dvar_sub,sub_tlef,vsize,'tf','K', &
            'Foliage canopy temperature','canopy_temperature',.true., &
            l_fill=.true.)
          sub_tlef_out => v2dvar_sub(sub_tlef)%rval
        end if
        if ( enable_sub2d_vars(sub_evp) ) then
          call setup_var(v2dvar_sub,sub_evp,vsize,'evspsbl','kg m-2 s-1', &
            'Total evapotranspiration flux','water_evaporation_flux',.true., &
            'time: mean')
          sub_evp_out => v2dvar_sub(sub_evp)%rval
        end if
        if ( enable_sub2d_vars(sub_scv) ) then
          call setup_var(v2dvar_sub,sub_scv,vsize,'snv','kg m-2', &
            'Liquid water equivalent of snow thickness', &
            'lwe_thickness_of_surface_snow_amount',.true.,'time: mean', &
            l_fill=.true.)
          sub_scv_out => v2dvar_sub(sub_scv)%rval
        end if
        if ( enable_sub2d_vars(sub_sena) ) then
          call setup_var(v2dvar_sub,sub_sena,vsize,'hfss','W m-2', &
            'Sensible heat flux','surface_upward_sensible_heat_flux', &
            .true.,'time: mean')
          sub_sena_out => v2dvar_sub(sub_sena)%rval
        end if

        vsize%k2 = 1
        v3dvar_sub(sub_u10m)%axis = 'xyw'
        v3dvar_sub(sub_v10m)%axis = 'xyw'
        v3dvar_sub(sub_t2m)%axis = 'xy2'
        v3dvar_sub(sub_q2m)%axis = 'xy2'
        if ( enable_sub3d_vars(sub_u10m) ) then
          call setup_var(v3dvar_sub,sub_u10m,vsize,'uas','m s-1', &
            '10 meter zonal wind component (westerly)', &
            'eastward_wind',.true.)
          sub_u10m_out => v3dvar_sub(sub_u10m)%rval
        end if
        if ( enable_sub3d_vars(sub_v10m) ) then
          call setup_var(v3dvar_sub,sub_v10m,vsize,'vas','m s-1', &
            '10 meter meridional wind component (southerly)', &
            'northward_wind',.true.)
          sub_v10m_out => v3dvar_sub(sub_v10m)%rval
        end if
        if ( enable_sub3d_vars(sub_t2m) ) then
          call setup_var(v3dvar_sub,sub_t2m,vsize,'tas','K', &
            '2 meter air temperature','air_temperature',.true.)
          sub_t2m_out => v3dvar_sub(sub_t2m)%rval
        end if
        if ( enable_sub3d_vars(sub_q2m) ) then
          call setup_var(v3dvar_sub,sub_q2m,vsize,'qas','1', &
            '2 meter air specific humidity','specific_humidity',.true.)
          sub_q2m_out => v3dvar_sub(sub_q2m)%rval
        end if
        vsize%k2 = 2
        v3dvar_sub(sub_smw)%axis = 'xys'
        v3dvar_sub(sub_runoff)%axis = 'xys'
        if ( enable_sub3d_vars(sub_smw) ) then
          call setup_var(v3dvar_sub,sub_smw,vsize,'mrso','kg m-2', &
            'Soil moisture content','soil_moisture_content', &
            .true.,l_fill=.true.)
          sub_smw_out => v3dvar_sub(sub_smw)%rval
        end if
        if ( enable_sub3d_vars(sub_runoff) ) then
          call setup_var(v3dvar_sub,sub_runoff,vsize,'mrro','kg m-2 day-1', &
            'Runoff flux','runoff_flux',.true.,'time: mean',l_fill=.true.)
          sub_runoff_out => v3dvar_sub(sub_runoff)%rval
        end if

        enable_sub_vars(1:nsub2dvars) = enable_sub2d_vars
        enable_sub_vars(nsub2dvars+1:nsubvars) = enable_sub3d_vars
        outstream(sub_stream)%nvar = countvars(enable_sub_vars,nsubvars)
        allocate(outstream(sub_stream)%ncvars%vlist(outstream(sub_stream)%nvar))
        outstream(sub_stream)%nfiles = 1
        allocate(outstream(sub_stream)%ncout(outstream(sub_stream)%nfiles))
        allocate(outstream(sub_stream)%cname_base(outstream(sub_stream)%nfiles))
        outstream(sub_stream)%cname_base(1) = 'SUB'

        vcount = 1
        do i = 1 , nsub2dvars
          if ( enable_sub_vars(i) ) then
            outstream(sub_stream)%ncvars%vlist(vcount)%vp => v2dvar_sub(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nsub3dvars
          if ( enable_sub_vars(i+nsub2dvars) ) then
            outstream(sub_stream)%ncvars%vlist(vcount)%vp => v3dvar_sub(i)
            vcount = vcount + 1
          end if
        end do
        outstream(sub_stream)%jl1 = vsize%j1
        outstream(sub_stream)%jl2 = vsize%j2
        outstream(sub_stream)%il1 = vsize%i1
        outstream(sub_stream)%il2 = vsize%i2
        outstream(sub_stream)%jg1 = joutsg1
        outstream(sub_stream)%jg2 = joutsg2
        outstream(sub_stream)%ig1 = ioutsg1
        outstream(sub_stream)%ig2 = ioutsg2
      end if

      if ( nstream == rad_stream ) then

        allocate(v2dvar_rad(nrad2dvars))
        allocate(v3dvar_rad(nrad3dvars))
        allocate(v4dvar_rad(nrad4dvars))
        enable_rad2d_vars = enable_rad_vars(1:nrad2dvars)
        enable_rad3d_vars = enable_rad_vars(nrad2dvars+1:nradvars-nrad4dvars)
        enable_rad4d_vars = enable_rad_vars(nrad2dvars+nrad3dvars+1:nradvars)
        ! This variables are always present

        call setup_common_vars(vsize,v2dvar_rad,rad_xlon, &
                  rad_xlat,rad_topo,rad_mask,rad_ps)

        ! The following may be enabled/disabled

        if ( enable_rad2d_vars(rad_frsa) ) then
          call setup_var(v2dvar_rad,rad_frsa,vsize,'rsns','W m-2', &
            'Surface net downward shortwave flux', &
            'surface_net_downward_shortwave_flux', .true.)
          rad_frsa_out => v2dvar_rad(rad_frsa)%rval
        end if
        if ( enable_rad2d_vars(rad_frla) ) then
          call setup_var(v2dvar_rad,rad_frla,vsize,'rsnl','W m-2', &
            'Surface net upward longwave flux', &
            'surface_net_upward_longwave_flux', .true.)
          rad_frla_out => v2dvar_rad(rad_frla)%rval
        end if
        if ( enable_rad2d_vars(rad_clrst) ) then
          call setup_var(v2dvar_rad,rad_clrst,vsize,'rtnscl','W m-2', &
            'Clearsky top of atmosphere net downward shortwave flux', &
            'toa_net_downward_shortwave_flux_assuming_clear_sky',.true.)
          rad_clrst_out => v2dvar_rad(rad_clrst)%rval
        end if
        if ( enable_rad2d_vars(rad_clrss) ) then
          call setup_var(v2dvar_rad,rad_clrss,vsize,'rsnscl','W m-2', &
            'Clearsky surface net downward shortwave flux', &
            'surface_net_downward_shortwave_flux_assuming_clear_sky',.true.)
          rad_clrss_out => v2dvar_rad(rad_clrss)%rval
        end if
        if ( enable_rad2d_vars(rad_clrlt) ) then
          call setup_var(v2dvar_rad,rad_clrlt,vsize,'rtnlcl','W m-2', &
            'Clearsky top of atmosphere net upward longwave flux', &
            'toa_net_upward_longwave_flux_assuming_clear_sky',.true.)
          rad_clrlt_out => v2dvar_rad(rad_clrlt)%rval
        end if
        if ( enable_rad2d_vars(rad_clrls) ) then
          call setup_var(v2dvar_rad,rad_clrls,vsize,'rsnlcl','W m-2', &
            'Clearsky net upward longwave flux', &
            'surface_net_upward_longwave_flux_assuming_clear_sky',.true.)
          rad_clrls_out => v2dvar_rad(rad_clrls)%rval
        end if
        if ( enable_rad2d_vars(rad_solin) ) then
          call setup_var(v2dvar_rad,rad_solin,vsize,'rts','W m-2', &
            'Top of atmosphere incoming shortwave flux', &
            'toa_incoming_shortwave_flux',.true.)
          rad_solin_out => v2dvar_rad(rad_solin)%rval
        end if
        if ( enable_rad2d_vars(rad_sabtp) ) then
          call setup_var(v2dvar_rad,rad_sabtp,vsize,'rsnt','W m-2', &
            'Net top of atmosphere upward shortwave flux', &
            'toa_net_upward_shortwave_flux',.true.)
          rad_sabtp_out => v2dvar_rad(rad_sabtp)%rval
        end if
        if ( enable_rad2d_vars(rad_totcf) ) then
          call setup_var(v2dvar_rad,rad_totcf,vsize,'clt','1', &
            'Total cloud fraction','cloud_area_fraction',.true.)
          rad_totcf_out => v2dvar_rad(rad_totcf)%rval
        end if
        if ( enable_rad2d_vars(rad_totcl) ) then
          call setup_var(v2dvar_rad,rad_totcl,vsize,'clwvi','kg m-2', &
            'Total columnar liquid water content', &
            'atmosphere_cloud_condensed_water_content',.true.)
          rad_totcl_out => v2dvar_rad(rad_totcl)%rval
        end if
        if ( enable_rad2d_vars(rad_totci) ) then
          call setup_var(v2dvar_rad,rad_totci,vsize,'clivi','kg m-2', &
            'Total columnar ice water content', &
            'atmosphere_ice_condensed_water_content',.true.)
          rad_totci_out => v2dvar_rad(rad_totci)%rval
        end if
        if ( enable_rad2d_vars(rad_firtp) ) then
          call setup_var(v2dvar_rad,rad_firtp,vsize,'rtl','W m-2', &
            'Top of atmosphere net upward longwave flux', &
            'toa_net_upward_longwave_flux',.true.)
          rad_firtp_out => v2dvar_rad(rad_firtp)%rval
        end if

        vsize%k2 = kz
        if ( enable_rad3d_vars(rad_cld) ) then
          call setup_var(v3dvar_rad,rad_cld,vsize,'cl','1', &
            'Cloud fractional cover', &
            'cloud_area_fraction_in_atmosphere_layer',.true.)
          rad_cld_out => v3dvar_rad(rad_cld)%rval
        end if
        if ( enable_rad3d_vars(rad_clwp) ) then
          call setup_var(v3dvar_rad,rad_clwp,vsize,'clw','g m-2', &
            'Cloud liquid water path','thickness_of_liquid_water_cloud',.true.)
          rad_clwp_out => v3dvar_rad(rad_clwp)%rval
        end if
        if ( enable_rad3d_vars(rad_qrs) ) then
          call setup_var(v3dvar_rad,rad_qrs,vsize,'qrs','K s-1', &
            'Shortwave radiation heating rate', &
            'tendency_of_air_temperature_due_to_shortwave_heating',.true.)
          rad_qrs_out => v3dvar_rad(rad_qrs)%rval
        end if
        if ( enable_rad3d_vars(rad_qrl) ) then
          call setup_var(v3dvar_rad,rad_qrl,vsize,'qrl','K s-1', &
            'Longwave radiation heating rate', &
            'tendency_of_air_temperature_due_to_longwave_heating',.true.)
          rad_qrl_out => v3dvar_rad(rad_qrl)%rval
        end if
        if ( idiag > 0 ) then
          if ( enable_rad4d_vars(rad_taucl) ) then
            v4dvar_rad(rad_taucl)%axis = 'xyzS'
            call setup_var(v4dvar_rad,rad_taucl,vsize,'taucl','1', &
              'Cloud liquid water optical depth', &
              'atmosphere_optical_thickness_due_to_cloud_liquid_water',.true.)
            rad_taucl_out => v4dvar_rad(rad_taucl)%rval
          end if
          if ( enable_rad4d_vars(rad_tauci) ) then
            v4dvar_rad(rad_tauci)%axis = 'xyzS'
            call setup_var(v4dvar_rad,rad_tauci,vsize,'tauci','1', &
              'Cloud ice optical depth', &
              'atmosphere_optical_thickness_due_to_cloud_ice',.true.)
            rad_tauci_out => v4dvar_rad(rad_tauci)%rval
          end if
        else
          enable_rad4d_vars(rad_taucl:rad_tauci) = .false.
        end if

        enable_rad_vars(1:nrad2dvars) = enable_rad2d_vars
        enable_rad_vars(nrad2dvars+1:nradvars-nrad4dvars) = enable_rad3d_vars
        enable_rad_vars(nrad2dvars+nrad3dvars+1:nradvars) = enable_rad4d_vars
        outstream(rad_stream)%nvar = countvars(enable_rad_vars,nradvars)
        allocate(outstream(rad_stream)%ncvars%vlist(outstream(rad_stream)%nvar))
        outstream(rad_stream)%nfiles = 1
        allocate(outstream(rad_stream)%ncout(outstream(rad_stream)%nfiles))
        allocate(outstream(rad_stream)%cname_base(outstream(rad_stream)%nfiles))
        outstream(rad_stream)%cname_base(1) = 'RAD'

        vcount = 1
        do i = 1 , nrad2dvars
          if ( enable_rad_vars(i) ) then
            outstream(rad_stream)%ncvars%vlist(vcount)%vp => v2dvar_rad(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nrad3dvars
          if ( enable_rad_vars(i+nrad2dvars) ) then
            outstream(rad_stream)%ncvars%vlist(vcount)%vp => v3dvar_rad(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nrad4dvars
          if ( enable_rad_vars(i+nrad2dvars+nrad3dvars) ) then
            outstream(rad_stream)%ncvars%vlist(vcount)%vp => v4dvar_rad(i)
            vcount = vcount + 1
          end if
        end do

        outstream(rad_stream)%opar%l_specint = .true.
        outstream(rad_stream)%jl1 = vsize%j1
        outstream(rad_stream)%jl2 = vsize%j2
        outstream(rad_stream)%il1 = vsize%i1
        outstream(rad_stream)%il2 = vsize%i2
        outstream(rad_stream)%jg1 = jout1
        outstream(rad_stream)%jg2 = jout2
        outstream(rad_stream)%ig1 = iout1
        outstream(rad_stream)%ig2 = iout2
      end if

      if ( nstream == lak_stream ) then

        allocate(v2dvar_lak(nlak2dvars))
        allocate(v3dvar_lak(nlak3dvars))
        enable_lak2d_vars = enable_lak_vars(1:nlak2dvars)
        enable_lak3d_vars = enable_lak_vars(nlak2dvars+1:nlakvars)

        ! This variables are always present

        call setup_common_vars(vsize,v2dvar_lak,lak_xlon, &
                  lak_xlat,lak_topo,lak_mask,lak_ps)

        ! The following may be enabled/disabled

        if ( enable_lak2d_vars(lak_tg) ) then
          call setup_var(v2dvar_lak,lak_tg,vsize,'ts','K', &
            'Ground temperature','surface_temperature',.true.)
          lak_tg_out => v2dvar_lak(lak_tg)%rval
        end if
        if ( enable_lak2d_vars(lak_tpr) ) then
          call setup_var(v2dvar_lak,lak_tpr,vsize,'pr','kg m-2 day-1', &
            'Total precipitation flux','precipitation_flux',.true.,'time: mean')
          lak_tpr_out => v2dvar_lak(lak_tpr)%rval
        end if
        if ( enable_lak2d_vars(lak_scv) ) then
          call setup_var(v2dvar_lak,lak_scv,vsize,'snv','kg m-2', &
            'Liquid water equivalent snow depth', &
            'lwe_thickness_of_surface_snow_amount',.true.,'time: mean', &
            l_fill=.true.)
          lak_scv_out => v2dvar_lak(lak_scv)%rval
        end if
        if ( enable_lak2d_vars(lak_sena) ) then
          call setup_var(v2dvar_lak,lak_sena,vsize,'hfss','W m-2', &
            'Sensible heat flux','surface_upward_sensible_heat_flux',.true.)
          lak_sena_out => v2dvar_lak(lak_sena)%rval
        end if
        if ( enable_lak2d_vars(lak_flw) ) then
          call setup_var(v2dvar_lak,lak_flw,vsize,'rsnl','W m-2', &
            'Net longwave energy flux','net_upward_longwave_flux_in_air',.true.)
          lak_flw_out => v2dvar_lak(lak_flw)%rval
        end if
        if ( enable_lak2d_vars(lak_fsw) ) then
          call setup_var(v2dvar_lak,lak_fsw,vsize,'rsns','W m-2', &
            'Net shortwave energy flux', 'net_downward_shortwave_flux_in_air', &
            .true.)
          lak_fsw_out => v2dvar_lak(lak_fsw)%rval
        end if
        if ( enable_lak2d_vars(lak_fld) ) then
          call setup_var(v2dvar_lak,lak_fld,vsize,'rsdl','W m-2', &
            'Downward longwave flux at surface in air', &
            'surface_downwelling_longwave_flux_in_air',.true.)
          lak_fld_out => v2dvar_lak(lak_fld)%rval
        end if
        if ( enable_lak2d_vars(lak_sina) ) then
          call setup_var(v2dvar_lak,lak_sina,vsize,'rsds','W m-2', &
            'Downward shortwave flux at surface in air', &
            'surface_downwelling_shortwave_flux_in_air',.true.)
          lak_sina_out => v2dvar_lak(lak_sina)%rval
        end if
        if ( enable_lak2d_vars(lak_aldirs) ) then
          call setup_var(v2dvar_lak,lak_aldirs,vsize,'aldirs','1', &
            'Surface albedo to direct shortwave radiation', &
            'surface_albedo_short_wave_direct',.true.)
          lak_aldirs_out => v2dvar_lak(lak_aldirs)%rval
        end if
        if ( enable_lak2d_vars(lak_aldifs) ) then
          call setup_var(v2dvar_lak,lak_aldifs,vsize,'aldifs','1', &
            'Surface albedo to diffuse shortwave radiation', &
            'surface_albedo_short_wave_diffuse',.true.)
          lak_aldifs_out => v2dvar_lak(lak_aldifs)%rval
        end if
        if ( enable_lak2d_vars(lak_evp) ) then
          call setup_var(v2dvar_lak,lak_evp,vsize,'evspsbl','mm s-1', &
            'Water evaporation flux', &
            'water_evaporation_flux_where_sea_ice',.true.)
          lak_evp_out => v2dvar_lak(lak_evp)%rval
        end if
        if ( enable_lak2d_vars(lak_ice) ) then
          call setup_var(v2dvar_lak,lak_ice,vsize,'lakice','mm', &
            'Lake Ice depth', 'seaice_depth',.true.)
          lak_ice_out => v2dvar_lak(lak_ice)%rval
        end if

        vsize%k2 = ndpmax
        v3dvar_lak(lak_tlake)%axis = 'xyd'
        if ( enable_lak3d_vars(lak_tlake) ) then
          call setup_var(v3dvar_lak,lak_tlake,vsize,'lakets','K', &
            'Lake water temperature','water_temperature',.true.,l_fill=.true.)
          lak_tlake_out => v3dvar_lak(lak_tlake)%rval
          lak_tlake_out = dmissval
        end if

        enable_lak_vars(1:nlak2dvars) = enable_lak2d_vars
        enable_lak_vars(nlak2dvars+1:nlakvars) = enable_lak3d_vars
        outstream(lak_stream)%nvar = countvars(enable_lak_vars,nlakvars)
        allocate(outstream(lak_stream)%ncvars%vlist(outstream(lak_stream)%nvar))
        outstream(lak_stream)%nfiles = 1
        allocate(outstream(lak_stream)%ncout(outstream(lak_stream)%nfiles))
        allocate(outstream(lak_stream)%cname_base(outstream(lak_stream)%nfiles))
        outstream(lak_stream)%cname_base(1) = 'LAK'

        vcount = 1
        do i = 1 , nlak2dvars
          if ( enable_lak_vars(i) ) then
            outstream(lak_stream)%ncvars%vlist(vcount)%vp => v2dvar_lak(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nlak3dvars
          if ( enable_lak_vars(i+nlak2dvars) ) then
            outstream(lak_stream)%ncvars%vlist(vcount)%vp => v3dvar_lak(i)
            vcount = vcount + 1
          end if
        end do
        outstream(lak_stream)%jl1 = vsize%j1
        outstream(lak_stream)%jl2 = vsize%j2
        outstream(lak_stream)%il1 = vsize%i1
        outstream(lak_stream)%il2 = vsize%i2
        outstream(lak_stream)%jg1 = jout1
        outstream(lak_stream)%jg2 = jout2
        outstream(lak_stream)%ig1 = iout1
        outstream(lak_stream)%ig2 = iout2
      end if

      if ( nstream == slaboc_stream ) then

        allocate(v2dvar_slaboc(nslaboc2dvars))
        allocate(v3dvar_slaboc(nslaboc3dvars))

        ! This variables are always present

        call setup_common_vars(vsize,v2dvar_slaboc,slab_xlon, &
                  slab_xlat,slab_topo,slab_mask,slab_ps)

        vsize%k2 = 12
        v3dvar_slaboc(slab_qflx)%axis = 'xyM'
        call setup_var(v3dvar_slaboc,slab_qflx,vsize,'qflx','W m-2', &
            'heat flux correction from slab ocean model', &
            'heat_flux_correction',.false.,'time: mean (interval 1 month)', &
            l_fill=.true.)
        slab_qflx_out => v3dvar_slaboc(slab_qflx)%rval

        outstream(slaboc_stream)%nvar = nslabocvars
        allocate( &
          outstream(slaboc_stream)%ncvars%vlist(outstream(slaboc_stream)%nvar))
        outstream(slaboc_stream)%nfiles = 1
        allocate( &
          outstream(slaboc_stream)%ncout(outstream(slaboc_stream)%nfiles))
        allocate( &
          outstream(slaboc_stream)%cname_base(outstream(slaboc_stream)%nfiles))
        outstream(slaboc_stream)%cname_base(1) = 'SOM'

        vcount = 1
        do i = 1 , nslaboc2dvars
          outstream(slaboc_stream)%ncvars%vlist(vcount)%vp => v2dvar_slaboc(i)
          vcount = vcount + 1
        end do
        do i = 1 , nslaboc3dvars
          outstream(slaboc_stream)%ncvars%vlist(vcount)%vp => v3dvar_slaboc(i)
          vcount = vcount + 1
        end do

        outstream(slaboc_stream)%jl1 = vsize%j1
        outstream(slaboc_stream)%jl2 = vsize%j2
        outstream(slaboc_stream)%il1 = vsize%i1
        outstream(slaboc_stream)%il2 = vsize%i2
        outstream(slaboc_stream)%jg1 = jout1
        outstream(slaboc_stream)%jg2 = jout2
        outstream(slaboc_stream)%ig1 = iout1
        outstream(slaboc_stream)%ig2 = iout2

      end if

      if ( nstream == opt_stream ) then

        allocate(v2dvar_opt(nopt2dvars))
        allocate(v3dvar_opt(nopt3dvars))
        enable_opt2d_vars = enable_opt_vars(1:nopt2dvars)
        enable_opt3d_vars = enable_opt_vars(nopt2dvars+1:noptvars)

        ! This variables are always present

        call setup_common_vars(vsize,v2dvar_opt,opt_xlon, &
                  opt_xlat,opt_topo,opt_mask,opt_ps)

        ! The following may be enabled/disabled

        if ( idirect > 0 ) then
          if ( enable_opt2d_vars(opt_acstoarf) ) then
            call setup_var(v2dvar_opt,opt_acstoarf,vsize,'acstoarf','W m-2', &
              'Top of atmosphere shortwave radiative forcing', &
              'toa_shortwave_radiative_forcing',.true.,'time: mean')
            opt_acstoarf_out => v2dvar_opt(opt_acstoarf)%rval
          end if
          if ( enable_opt2d_vars(opt_acstsrrf) ) then
            call setup_var(v2dvar_opt,opt_acstsrrf,vsize,'acstsrrf','W m-2', &
              'Surface shortwave radiative forcing', &
              'surface_shortwave_radiative_forcing',.true.,'time: mean')
            opt_acstsrrf_out => v2dvar_opt(opt_acstsrrf)%rval
          end if
          if ( enable_opt2d_vars(opt_acstalrf) ) then
            call setup_var(v2dvar_opt,opt_acstalrf,vsize,'acstalrf','W m-2', &
              'Top of atmosphere longwave radiative forcing', &
              'toa_longwave_radiative_forcing',.true.,'time: mean')
            opt_acstalrf_out => v2dvar_opt(opt_acstalrf)%rval
          end if
          if ( enable_opt2d_vars(opt_acssrlrf) ) then
            call setup_var(v2dvar_opt,opt_acssrlrf,vsize,'acssrlrf','W m-2', &
              'Surface longwave radiative forcing' , &
              'surface_longwave_radiative_forcing',.true.,'time: mean')
            opt_acssrlrf_out => v2dvar_opt(opt_acssrlrf)%rval
          end if
        else
          enable_opt2d_vars(opt_acstoarf:opt_acssrlrf) = .false.
        end if
        if ( enable_opt2d_vars(opt_aod) ) then
          call setup_var(v2dvar_opt,opt_aod,vsize,'aod','1', &
            'Aerosol optical thickness in the visible band' , &
            'atmosphere_optical_thickness_due_to_aerosol',.true.)
          opt_aod_out => v2dvar_opt(opt_aod)%rval
        end if

        vsize%k2 = kz
        if ( enable_opt3d_vars(opt_aext8) ) then
          call setup_var(v3dvar_opt,opt_aext8,vsize,'aext8','1', &
            'Aerosol optical depth', &
            'atmosphere_optical_thickness_due_to_aerosol',.true.)
          opt_aext8_out => v3dvar_opt(opt_aext8)%rval
        end if
        if ( enable_opt3d_vars(opt_assa8) ) then
          call setup_var(v3dvar_opt,opt_assa8,vsize,'assa8','1', &
            'Aerosol single scattering albedo', &
            'aerosol_single_scattering_albedo',.true.)
          opt_assa8_out => v3dvar_opt(opt_assa8)%rval
        end if
        if ( enable_opt3d_vars(opt_agfu8) ) then
          call setup_var(v3dvar_opt,opt_agfu8,vsize,'agfu8','1', &
            'Aerosol asymmetry parameter', &
            'aerosol_asymmetry_parameter',.true.)
          opt_agfu8_out => v3dvar_opt(opt_agfu8)%rval
        end if

        enable_opt_vars(1:nopt2dvars) = enable_opt2d_vars
        enable_opt_vars(nopt2dvars+1:noptvars) = enable_opt3d_vars
        outstream(opt_stream)%nvar = countvars(enable_opt_vars,noptvars)
        allocate(outstream(opt_stream)%ncvars%vlist(outstream(opt_stream)%nvar))
        outstream(opt_stream)%nfiles = 1
        allocate(outstream(opt_stream)%ncout(outstream(opt_stream)%nfiles))
        allocate(outstream(opt_stream)%cname_base(outstream(opt_stream)%nfiles))
        outstream(opt_stream)%cname_base(1) = 'OPT'

        vcount = 1
        do i = 1 , nopt2dvars
          if ( enable_opt_vars(i) ) then
            outstream(opt_stream)%ncvars%vlist(vcount)%vp => v2dvar_opt(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nopt3dvars
          if ( enable_opt_vars(i+nopt2dvars) ) then
            outstream(opt_stream)%ncvars%vlist(vcount)%vp => v3dvar_opt(i)
            vcount = vcount + 1
          end if
        end do
        outstream(opt_stream)%jl1 = vsize%j1
        outstream(opt_stream)%jl2 = vsize%j2
        outstream(opt_stream)%il1 = vsize%i1
        outstream(opt_stream)%il2 = vsize%i2
        outstream(opt_stream)%jg1 = jout1
        outstream(opt_stream)%jg2 = jout2
        outstream(opt_stream)%ig1 = iout1
        outstream(opt_stream)%ig2 = iout2

      end if

      if ( nstream == che_stream ) then

        allocate(v2dvar_che(nche2dvars))
        allocate(v3dvar_che(nche3dvars))
        enable_che2d_vars = enable_che_vars(1:nche2dvars)
        enable_che3d_vars = enable_che_vars(nche2dvars+1:nchevars)

        ! This variables are always present

        call setup_common_vars(vsize,v2dvar_che,che_xlon, &
                  che_xlat,che_topo,che_mask,che_ps)

        ! The following may be enabled/disabled

        if ( enable_che2d_vars(che_wdrflx) ) then
          call setup_var(v2dvar_che,che_wdrflx,vsize,'wdrflx', &
            'mg m-2 day-1','Wet deposition flux due to rainout', &
            'wet_deposition_flux_from_rainout',.true.,'time: mean')
          che_wdrflx_out => v2dvar_che(che_wdrflx)%rval
        end if
        if ( enable_che2d_vars(che_wdcflx) ) then
          call setup_var(v2dvar_che,che_wdcflx,vsize,'wdwflx', &
            'mg m-2 day-1','Wet deposition flux due to washout', &
            'wet_deposition_flux_from_washout',.true.,'time: mean')
          che_wdcflx_out => v2dvar_che(che_wdcflx)%rval
        end if
        if ( enable_che2d_vars(che_ddflx) ) then
          call setup_var(v2dvar_che,che_ddflx,vsize,'ddflx', &
            'mg m-2 day-1','Dry deposition flux', &
            'dry_deposition_flux',.true.,'time: mean')
          che_ddflx_out => v2dvar_che(che_ddflx)%rval
        end if
        if ( enable_che2d_vars(che_emflx) ) then
          call setup_var(v2dvar_che,che_emflx,vsize,'emflx', &
            'mg m-2 day-1','Surface emission flux', &
            'surface_emission_flux',.true.,'time: mean')
          che_emflx_out => v2dvar_che(che_emflx)%rval
        end if
        if ( enable_che2d_vars(che_ddvel) ) then
          call setup_var(v2dvar_che,che_ddvel,vsize,'ddvel', &
            'm s-1','Dry deposition velocity', &
            'dry_deposition_velocity',.true.,'time: mean')
          che_ddvel_out => v2dvar_che(che_ddvel)%rval
        end if
        if ( enable_che2d_vars(che_burden) ) then
          call setup_var(v2dvar_che,che_burden,vsize,'burden', &
            'mg m-2','Tracer total burden', &
            'tracer_burden',.true.,'time: mean')
          che_burden_out => v2dvar_che(che_burden)%rval
        end if

        vsize%k2 = kz
        if ( enable_che3d_vars(che_mixrat) ) then
          call setup_var(v3dvar_che,che_mixrat,vsize,'mixrat','1', &
            'Atmosphere tracer mixing ratio', &
            'atmosphere_mixing_ratio_of_tracer',.true.)
          che_mixrat_out => v3dvar_che(che_mixrat)%rval
        end if

        if ( ichdiag == 1 ) then
          if ( enable_che3d_vars(che_cheten) ) then
            call setup_var(v3dvar_che,che_cheten,vsize,'cheten', &
              'kg kg-1 s-1', 'Tendency of tracer due to chemical reactions', &
              'tendency_of_mixing_ratio_due_to_chemical_prod_loss',.true.)
            che_cheten_out => v3dvar_che(che_cheten)%rval
          end if
          if ( enable_che3d_vars(che_advhten) ) then
            call setup_var(v3dvar_che,che_advhten,vsize,'advhten', &
              'kg kg-1 s-1', 'Tendency of tracer due to horizontal advection', &
              'tendency_of_mixing_ratio_due_to_horizontal_advection',.true.)
            che_advhten_out => v3dvar_che(che_advhten)%rval
          end if
          if ( enable_che3d_vars(che_advvten) ) then
            call setup_var(v3dvar_che,che_advvten,vsize,'advvten', &
              'kg kg-1 s-1', 'Tendency of tracer due to vertical advection', &
              'tendency_of_mixing_ratio_due_to_vertical_advection',.true.)
            che_advvten_out => v3dvar_che(che_advvten)%rval
          end if
          if ( enable_che3d_vars(che_difhten) ) then
            call setup_var(v3dvar_che,che_difhten,vsize,'difhten', &
              'kg kg-1 s-1', 'Tendency of tracer due to diffusion', &
              'tendency_of_mixing_ratio_due_to_diffusion',.true.)
            che_difhten_out => v3dvar_che(che_difhten)%rval
          end if
          if ( enable_che3d_vars(che_cuten) ) then
            call setup_var(v3dvar_che,che_cuten,vsize,'cuten', &
              'kg kg-1 s-1', 'Tendency of tracer due to convective transport', &
              'tendency_of_mixing_ratio_due_to_convective_transport',.true.)
            che_cuten_out => v3dvar_che(che_cuten)%rval
          end if
          if ( enable_che3d_vars(che_tuten) ) then
            call setup_var(v3dvar_che,che_tuten,vsize,'tuten', &
              'kg kg-1 s-1', 'Tendency of tracer due to vertical turbolence', &
              'tendency_of_mixing_ratio_due_to_vertical_turbolence',.true.)
            che_tuten_out => v3dvar_che(che_tuten)%rval
          end if
          if ( enable_che3d_vars(che_raiten) ) then
            call setup_var(v3dvar_che,che_raiten,vsize,'raiten', &
              'kg kg-1 s-1', 'Tendency of tracer due to total rainout', &
              'tendency_of_mixing_ratio_due_to_total_rainout',.true.)
            che_raiten_out => v3dvar_che(che_raiten)%rval
          end if
          if ( enable_che3d_vars(che_wasten) ) then
            call setup_var(v3dvar_che,che_wasten,vsize,'wasten', &
              'kg kg-1 s-1', 'Tendency of tracer due to total washout', &
              'tendency_of_mixing_ratio_due_to_total_washout',.true.)
            che_wasten_out => v3dvar_che(che_wasten)%rval
          end if
          if ( enable_che3d_vars(che_bdyten) ) then
            call setup_var(v3dvar_che,che_bdyten,vsize,'bdyten', &
              'kg kg-1 s-1', 'Tendency of tracer due to boundary conditions', &
              'tendency_of_mixing_ratio_due_to_boundary_conditions',.true.)
            che_bdyten_out => v3dvar_che(che_bdyten)%rval
          end if
          if ( enable_che3d_vars(che_sedten) ) then
            call setup_var(v3dvar_che,che_sedten,vsize,'sedten', &
              'kg kg-1 s-1', 'Tendency of tracer due to sedimentation', &
              'tendency_of_mixing_ratio_due_to_sedimentation',.true.)
            che_sedten_out => v3dvar_che(che_sedten)%rval
          end if
          if ( enable_che3d_vars(che_emten) ) then
            call setup_var(v3dvar_che,che_emten,vsize,'emiten', &
              'kg kg-1 s-1', 'Tendency of tracer due to emission', &
              'tendency_of_mixing_ratio_due_to_emission',.true.)
            che_emten_out => v3dvar_che(che_emten)%rval
          end if
          if ( ibltyp == 2 .or. ibltyp == 99 ) then
            if ( enable_che2d_vars(che_pblten) ) then
              call setup_var(v2dvar_che,che_pblten,vsize,'pblten', &
                'kg kg-1 s-1', 'Tendency of tracer due to UW PBL', &
                'tendency_of_mixing_ratio_due_to_uw_pbl',.true.)
              che_pblten_out => v2dvar_che(che_pblten)%rval
            end if
          else
            enable_che2d_vars(che_pblten) = .false.
          end if
        else
          enable_che3d_vars(che_cheten:che_emten) = .false.
          enable_che2d_vars(che_pblten) = .false.
        end if
        enable_che_vars(1:nche2dvars) = enable_che2d_vars
        enable_che_vars(nche2dvars+1:nchevars) = enable_che3d_vars

        outstream(che_stream)%nvar = countvars(enable_che_vars,nchevars)
        allocate(outstream(che_stream)%ncvars%vlist(outstream(che_stream)%nvar))
        outstream(che_stream)%nfiles = ntr
        allocate(outstream(che_stream)%ncout(outstream(che_stream)%nfiles))
        allocate(outstream(che_stream)%cname_base(outstream(che_stream)%nfiles))
        do itr = 1 , ntr
          outstream(che_stream)%cname_base(itr) = chtrname(itr)
        end do

        vcount = 1
        do i = 1 , nche2dvars
          if ( enable_che_vars(i) ) then
            outstream(che_stream)%ncvars%vlist(vcount)%vp => v2dvar_che(i)
            vcount = vcount + 1
          end if
        end do
        do i = 1 , nche3dvars
          if ( enable_che_vars(i+nche2dvars) ) then
            outstream(che_stream)%ncvars%vlist(vcount)%vp => v3dvar_che(i)
            vcount = vcount + 1
          end if
        end do
        outstream(che_stream)%jl1 = vsize%j1
        outstream(che_stream)%jl2 = vsize%j2
        outstream(che_stream)%il1 = vsize%i1
        outstream(che_stream)%il2 = vsize%i2
        outstream(che_stream)%jg1 = jout1
        outstream(che_stream)%jg2 = jout2
        outstream(che_stream)%ig1 = iout1
        outstream(che_stream)%ig2 = iout2

      end if

      outstream(nstream)%opar%pname = 'RegCM Model'
      outstream(nstream)%opar%l_sync = lsync
      outstream(nstream)%opar%l_band = (i_band == 1)

      if ( parallel_out ) then
        outstream(nstream)%opar%mpi_comm = mycomm
        outstream(nstream)%opar%mpi_info = ncout_mpi_info
#ifdef NETCDF4_HDF5
        outstream(nstream)%opar%mpi_iotype = nf90_mpiio
#endif
        ! The "global" indexes in the output stream refer to the INTERNAL
        ! CROSS grid, i.e. for processor 0 this is (2,2) => (1,1) so we must
        ! subtract 1 line/column to rebase on pixel (2,2) of the internal model
        ! cross points grid to point (1,1).
        outstream(nstream)%opar%global_jstart = global_out_jstart
        outstream(nstream)%opar%global_jend   = global_out_jend
        outstream(nstream)%opar%global_istart = global_out_istart
        outstream(nstream)%opar%global_iend   = global_out_iend
        if ( nstream == sub_stream ) then
          outstream(nstream)%opar%global_jstart = &
            (outstream(nstream)%opar%global_jstart-1)*nsg+1
          outstream(nstream)%opar%global_jend =   &
            outstream(nstream)%opar%global_jend*nsg
          outstream(nstream)%opar%global_istart = &
            (outstream(nstream)%opar%global_istart-1)*nsg+1
          outstream(nstream)%opar%global_iend =   &
            outstream(nstream)%opar%global_iend*nsg
        end if
      end if

    end do enabled_stream_loop

    ! Allocate space to collect all from all CPUs if not parallel output

    if ( .not. parallel_out ) then
      if ( myid == iocpu ) then
        kkz = 2
        n4dd = 4
        if ( atm_stream > 0 .or. rad_stream > 0 ) then
          kkz = max(kz,kkz)
        end if
        if ( lak_stream > 0 ) then
          kkz = max(ndpmax,kkz)
        end if
        call getmem2d(io2d,jout1,jout2,iout1,iout2,'ncout:io2d')
        call getmem3d(io3d,jout1,jout2,iout1,iout2,1,kkz,'ncout:io3d')
        call getmem4d(io4d,jout1,jout2,iout1,iout2,1,kkz,1,n4dd,'ncout:io4d')
        if ( sub_stream > 0 ) then
          call getmem2d(io2dsg,joutsg1,joutsg2,ioutsg1,ioutsg2,'ncout:io2dsg')
          call getmem3d(io3dsg,joutsg1,joutsg2,ioutsg1,ioutsg2, &
                        1,2,'ncout:io3dsg')
        end if
      end if
    end if

  end subroutine init_output_streams

  integer(ik4) function countvars(eflags,ntot)
    implicit none
    integer(ik4) , intent(in) :: ntot
    logical , dimension(ntot) , intent(in) :: eflags
    integer(ik4) :: i
    countvars = nbase
    do i = nbase+1 , ntot
      if ( eflags(i) ) countvars = countvars + 1
    end do
  end function countvars

  subroutine newoutfiles(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=36) :: fbname
    character(len=36) :: cdate
    class(ncvariable_standard) , pointer :: vp
    real(rk8) , pointer , dimension(:,:) :: tmp2d
    real(rk8) , pointer , dimension(:,:) :: pnt2d => null()
    integer(ik4) :: i , j , ivar
    real(rk8) :: dummy

    if ( .not. parallel_out .and. myid /= iocpu ) then
      stream_loop: &
      do i = 1 , maxstreams
        file_loop: &
        do j = 1 , outstream(i)%nfiles
          var_loop: &
          do ivar = 1 , outstream(i)%nvar
            vp => outstream(i)%ncvars%vlist(ivar)%vp
            select type(vp)
              type is (ncvariable2d_real)
                if ( vp%lrecords ) cycle var_loop
                call grid_collect(vp%rval,pnt2d, &
                                  vp%j1,vp%j2,vp%i1,vp%i2,outstream(i)%l_sub)
              class default
                cycle var_loop
            end select
          end do var_loop
        end do file_loop
      end do stream_loop
      return
    end if

    stream_loop_par: &
    do i = 1 , maxstreams

      file_loop_par: &
      do j = 1 , outstream(i)%nfiles

        if ( i == slaboc_stream ) then
          write (fbname,'(a,a,a)') trim(outstream(i)%cname_base(j)) , &
            '.YYYYMMDDHH'
          outstream(i)%opar%fname = &
            trim(dirglob)//pthsep//trim(domname)//'_'//trim(fbname)//'.nc'
          outstream(i)%opar%zero_date = idate
        else
          write (fbname,'(a,a,i10)') trim(outstream(i)%cname_base(j)) , &
            '.', toint10(idate)
          outstream(i)%opar%fname = &
            trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)//'.nc'
          outstream(i)%opar%zero_date = idate
        end if

        if ( myid == italk ) then
          write(stdout,*) 'Opening new output file ', &
            trim(outstream(i)%opar%fname)
        end if
        call outstream_setup(outstream(i)%ncout(j),outstream(i)%opar)

        ! Initial and Boundary data sources

        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_string('model_icbc_data_source',dattyp))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_string('model_sst_data_source',ssttyp))

        ! Buffer Zone Control relaxation + diffusion term params

        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('boundary_nspgx',nspgx))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('boundary_nspgd',nspgd))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_real8('boundary_high_nudge',high_nudge))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_real8('boundary_medium_nudge',medium_nudge))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_real8('boundary_low_nudge',low_nudge))

        ! Perturbation control for ensembles

        if ( ensemble_run ) then
          if ( lperturb_topo ) then
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('perturbation_topo_percent',perturb_frac_topo))
          end if
          if ( lperturb_ts ) then
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('perturbation_ts_percent',perturb_frac_ts))
          end if
          if ( lperturb_ps ) then
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('perturbation_ps_percent',perturb_frac_ps))
          end if
          if ( lperturb_t ) then
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('perturbation_t_percent',perturb_frac_t))
          end if
          if ( lperturb_u ) then
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('perturbation_u_percent',perturb_frac_u))
          end if
          if ( lperturb_v ) then
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('perturbation_v_percent',perturb_frac_v))
          end if
          if ( lperturb_q ) then
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('perturbation_q_percent',perturb_frac_q))
          end if
        end if

        ! Model start/restart control

        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_logical('model_is_restarted',ifrest))
        cdate = tochar(idate0)
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_string('model_simulation_initial_start',cdate))
        cdate = tochar(idate1)
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_string('model_simulation_start',cdate))
        cdate = tochar(idate2)
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_string('model_simulation_end',cdate))

        ! Model timing parameters

        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_real8('atmosphere_time_step_in_seconds',dt))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_real8('surface_interaction_time_step_in_seconds',dtsrf))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_real8('radiation_scheme_time_step_in_minuts',dtrad))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_real8('absorption_emission_time_step_in_hours',dtabem))

        ! Model Physics

        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('lateral_boundary_condition_scheme',iboudy))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('semi_lagrangian_advection_scheme',isladvec))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('boundary_layer_scheme',ibltyp))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('cumulus_convection_scheme_lnd',icup_lnd))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('cumulus_convection_scheme_ocn',icup_ocn))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('moisture_scheme',ipptls))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('ocean_flux_scheme',iocnflx))
        if ( iocnflx == 2 ) then
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('zeng_ocean_roughness_formula',iocnrough))
        end if
        if ( iocncpl == 1 ) then
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('coupled_ocean_run',iocncpl))
        end if
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('pressure_gradient_scheme',ipgf))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('surface_emissivity_factor_computed',iemiss))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('lake_model_activated',lakemod))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('chemical_aerosol_scheme_activated',ichem))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_string('ipcc_scenario_code',scenario))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('diurnal_cycle_sst_scheme',idcsst))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('simple_sea_ice_scheme',iseaice))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('seasonal_desert_albedo',idesseas))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('convective_lwp_as_large_scale',iconvlwp))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('rrtm_radiation_scheme_activated',irrtm))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('climatic_ozone_input_dataset',iclimao3))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('static_solar_constant_used',isolconst))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('cumulus_cloud_model',icumcloud))
        if ( ipptls == 1 ) then
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('subex_bottom_level_with_no_clouds',ncld))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_auto_conversion_rate_for_land',qck1land))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_auto_conversion_rate_for_ocean',qck1oce))
          call outstream_addatt(outstream(i)%ncout(j), &
           ncattribute_real8('subex_gultepe_factor_when_rain_for_land',gulland))
          call outstream_addatt(outstream(i)%ncout(j), &
           ncattribute_real8('subex_gultepe_factor_when_rain_for_ocean',guloce))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_rh_with_fcc_one',rhmax))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_rh_threshold_for_land',rh0land))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_rh_threshold_for_ocean',rh0oce))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_limit_temperature',tc0))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_land_raindrop_evaporation_rate',cevaplnd))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_ocean_raindrop_evaporation_rate',cevapoce))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_land_raindrop_accretion_rate',caccrlnd))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_ocean_raindrop_accretion_rate',caccroce))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_cloud_fraction_maximum',cftotmax))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_condensation_threshold',conf))
          call outstream_addatt(outstream(i)%ncout(j), &
         ncattribute_real8('subex_cloud_fraction_max_for_convection',clfrcvmax))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('subex_cloud_liqwat_max_for_convection',cllwcv))
        else if ( ipptls == 2 ) then
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_logical('micro_statistics',stats))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_logical('micro_budget_verification',budget_compute))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('micro_super_saturation_option',nssopt))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('micro_autoconversion_option',kautoconv))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('micro_semi_implicit_option',ksemi))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('micro_fall_speed_rain',vqxr))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('micro_fall_speed_ice',vqxi))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('micro_fall_speed_snow',vqxs))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('micro_autoconv_khair',zauto_rate_khair))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('micro_autoconv_kessl',zauto_rate_kessl))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('micro_autoconv_klepi',zauto_rate_klepi))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('micro_autoconv_timescale_sund',rkconv))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('micro_min_cloud_coverage',rcovpmin))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('micro_evap_rate',rpecons))
        end if
        if ( any(icup == 2) ) then
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('grell_scheme_closure',igcc))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_min_shear_on_precip',shrmin))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_max_shear_on_precip',shrmax))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_min_precip_efficiency',edtmin))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_max_precip_efficiency',edtmax))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_min_precip_efficiency_o',edtmino))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_max_precip_efficiency_o',edtmaxo))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_min_precip_efficiency_x',edtminx))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_max_precip_efficiency_x',edtmaxx))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_min_shear_on_precip_on_ocean',shrmin_ocn))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_max_shear_on_precip_on_ocean',shrmax_ocn))
          call outstream_addatt(outstream(i)%ncout(j), &
           ncattribute_real8('grell_min_precip_efficiency_on_ocean',edtmin_ocn))
          call outstream_addatt(outstream(i)%ncout(j), &
           ncattribute_real8('grell_max_precip_efficiency_on_ocean',edtmax_ocn))
          call outstream_addatt(outstream(i)%ncout(j), &
        ncattribute_real8('grell_min_precip_efficiency_o_on_ocean',edtmino_ocn))
          call outstream_addatt(outstream(i)%ncout(j), &
        ncattribute_real8('grell_max_precip_efficiency_o_on_ocean',edtmaxo_ocn))
          call outstream_addatt(outstream(i)%ncout(j), &
        ncattribute_real8('grell_min_precip_efficiency_x_on_ocean',edtminx_ocn))
          call outstream_addatt(outstream(i)%ncout(j), &
        ncattribute_real8('grell_max_precip_efficiency_x_on_ocean',edtmaxx_ocn))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_max_depth_of_stable_layer',pbcmax))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_min_depth_of_cloud',mincld))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_min_convective_heating',htmin))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_max_convective_heating',htmax))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_max_cloud_base_height',skbmax))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('grell_FC_ABE_removal_timescale',dtauc))
         end if
        if ( any(icup == 4) ) then
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_lowest_convection_sigma',minsig))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8( &
            'mit_autoconversion_threshold_mixing_over_ocean',elcrit_ocn))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8( &
            'mit_autoconversion_threshold_mixing_over_land',elcrit_lnd))
          call outstream_addatt(outstream(i)%ncout(j), &
           ncattribute_real8('mit_autoconversion_threshold_temperature',tlcrit))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_mixing_coefficient_in_entrainment',entp))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_fractional_area_unsaturated_downdraft',sigd))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_fractional_precip_outside_cloud',sigs))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_pressure_velocity_of_rain',omtrain))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_pressure_velocity_of_snow',omtsnow))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_rain_evaporation_coefficient',coeffr))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_snow_evaporation_coefficient',coeffs))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_momentum_transport_coefficient',cu))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_downdraft_velocity_coefficient',betae))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_max_parcel_neg_temp_perturbation',dtmax))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_approach_rate_quasi_eq_coeff_a',alphae))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_approach_rate_quasi_eq_coeff_d',damp))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_maximum_land_precipitation_efficiency', &
            epmax_lnd))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('mit_maximum_ocean_precipitation_efficiency', &
            epmax_ocn))
        end if
        if ( any(icup == 5) ) then
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('tiedtke_actual_scheme',iconv))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('tiedtke_entrainment_rate_downdraft',entrdd))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('tiedtke_entrainment_rate_deep',entrpen))
          if ( iconv == 4 ) then
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_detrainment_rate_deep',detrpen))
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_shallow_entrainment',entshalp))
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_cloud_cover_evap_over_land',rcuc_lnd))
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_cloud_cover_evap_over_ocean',rcuc_ocn))
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_coeff_evap_over_land',rcpec_lnd))
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_coeff_evap_over_ocean',rcpec_ocn))
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_critical_rh_over_land',rhebc_lnd))
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_critical_rh_over_ocean',rhebc_ocn))
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_cloud_water_conv_over_land',rprc_lnd))
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_cloud_water_conv_over_ocean',rprc_ocn))
          else
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_max_entrainment',entrmax))
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_entrainment_rate_shallow',entrscv))
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_entrainment_rate_midlevel',entrmid))
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tiedtke_conversion_coefficient',cprcon))
          end if
        end if
        if ( any(icup == 6) ) then
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('kf_trigger',kf_trigger))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('kf_entrainment_rate',kf_entrate))
        end if
        if ( ibltyp == 1 .or. ibltyp == 99 ) then
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('holtslag_critical_ocean_richardson',ricr_ocn))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('holtslag_critical_land_richardson',ricr_lnd))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('holtslag_zhnew_factor',zhnew_fac))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('holtslag_th10_estimate',ifaholtth10))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('holtslag_th10_maximize',ifaholtmax))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('holtslag_th10_minimize',ifaholtmin))
        end if
        if ( ibltyp == 2 .or. ibltyp == 99 ) then
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('uwpbl_advection_scheme',iuwvadv))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('uwpbl_cloud_evap_entr_incr_efficiency',atwo))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('uwpbl_eddy_LS_stable_PBL_scaling',rstbl))
        end if
        if ( irrtm == 1 ) then
          call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('rrtm_opt_properties_calculation_sw_vap',inflgsw))
          call outstream_addatt(outstream(i)%ncout(j), &
         ncattribute_integer('rrtm_opt_properties_calculation_sw_ice',iceflgsw))
          call outstream_addatt(outstream(i)%ncout(j), &
         ncattribute_integer('rrtm_opt_properties_calculation_sw_liq',iceflgsw))
          call outstream_addatt(outstream(i)%ncout(j), &
         ncattribute_integer('rrtm_opt_properties_calculation_lw_vap',inflglw))
          call outstream_addatt(outstream(i)%ncout(j), &
         ncattribute_integer('rrtm_opt_properties_calculation_lw_ice',iceflglw))
          call outstream_addatt(outstream(i)%ncout(j), &
         ncattribute_integer('rrtm_opt_properties_calculation_lw_liq',iceflglw))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('rrtm_idrv',idrv))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('rrtm_cloud_overlap_hypothesis',icld))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('rrtm_random_number_generator',irng))
        end if
        if ( ichem == 1 ) then
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_string('chem_simulation_type',chemsimtype))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('chem_activate_reaction_solver',ichsolver))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('chem_activate_emission',ichsursrc))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('chem_activate_dry_deposition',ichdrdepo))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('chem_activate_convective_transport',ichcumtra))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('chem_activate_wet_removal',ichremlsc))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_integer('chem_dust_PSD_scheme',ichdustemd))
          call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('chem_enable_aerosol_radiation_feedback',idirect))
          call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('chem_enable_sulfate_indirect_effect',iindirect))
        end if
#ifdef CLM
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('clm_land_surface_dataset_selection', imask))
        call outstream_addatt(outstream(i)%ncout(j), &
          ncattribute_integer('clm_use_modified_lawrence_albedo', &
              ilawrence_albedo))
#endif
        if ( iocncpl == 1 ) then
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('cpl_coupling_timestep_in_seconds',cpldt))
        end if

        if ( islab_ocean == 1 ) then
          if ( do_qflux_adj ) &
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_logical('slabocean_qflux_adjusted_run',do_qflux_adj))
          if ( do_restore_sst ) &
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_logical('slabocean_do_restore_sst',do_restore_sst))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('slabocean_sst_restore_timescale', &
              sst_restore_timescale))
          call outstream_addatt(outstream(i)%ncout(j), &
            ncattribute_real8('slabocean_mixed_layer_depth', &
              mixed_layer_depth))
        end if

        if ( itweak == 1 ) then
          if ( itweak_sst == 1 ) then
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tweak_sst',sst_tweak))
          end if
          if ( itweak_temperature == 1 ) then
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tweak_temperature',temperature_tweak))
          end if
          if ( itweak_solar_irradiance == 1 ) then
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tweak_solar_irradiance',solar_tweak))
          end if
          if ( itweak_greenhouse_gases == 1 ) then
            dummy = gas_tweak_factors(1)
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tweak_co2_factor',dummy))
            dummy = gas_tweak_factors(2)
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tweak_ch4_factor',dummy))
            dummy = gas_tweak_factors(3)
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tweak_n2o_factor',dummy))
            dummy = gas_tweak_factors(4)
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tweak_cfc11_factor',dummy))
            dummy = gas_tweak_factors(5)
            call outstream_addatt(outstream(i)%ncout(j), &
              ncattribute_real8('tweak_cfc12_factor',dummy))
          end if
        end if

        do ivar = 1 , outstream(i)%nvar
          vp => outstream(i)%ncvars%vlist(ivar)%vp
          call outstream_addvar(outstream(i)%ncout(j),vp)
        end do

        call outstream_enable(outstream(i)%ncout(j),hsigma)

        if ( .not. parallel_out ) then
          if ( outstream(i)%l_sub ) then
            pnt2d => io2dsg
          else
            pnt2d => io2d
          end if
        end if

        var_loop_par: &
        do ivar = 1 , outstream(i)%nvar
          vp => outstream(i)%ncvars%vlist(ivar)%vp

          ! Collect on temp storage
          if ( .not. parallel_out ) then
            select type(vp)
              type is (ncvariable2d_real)
                if ( vp%lrecords ) cycle var_loop_par
                call grid_collect(vp%rval,pnt2d, &
                                  vp%j1,vp%j2,vp%i1,vp%i2,outstream(i)%l_sub)
                vp%j1 = outstream(i)%jg1
                vp%j2 = outstream(i)%jg2
                vp%i1 = outstream(i)%ig1
                vp%i2 = outstream(i)%ig2
                tmp2d => vp%rval
                vp%rval => pnt2d
              class default
                cycle var_loop_par
            end select
          else
            select type(vp)
              type is (ncvariable2d_real)
                if ( vp%lrecords ) cycle var_loop_par
              class default
                cycle var_loop_par
            end select
          end if

          call outstream_writevar(outstream(i)%ncout(j),vp)

          ! Reset pointer and shape
          if ( .not. parallel_out ) then
            select type(vp)
              type is (ncvariable2d_real)
                vp%rval => tmp2d
                vp%j1 = outstream(i)%jl1
                vp%j2 = outstream(i)%jl2
                vp%i1 = outstream(i)%il1
                vp%i2 = outstream(i)%il2
              class default
                cycle var_loop_par
            end select
          end if

        end do var_loop_par
      end do file_loop_par
    end do stream_loop_par
  end subroutine newoutfiles

  subroutine setup_common_vars(vsize,var,xlon,xlat,topo,mask,ps)
    implicit none
    type(varspan) , intent(in) :: vsize
    type(ncvariable2d_real) , dimension(:) , intent(inout) :: var
    integer(ik4), intent(in) :: xlon , xlat , topo , mask , ps
    if ( associated(xlon_out) ) then
      call setup_var(var,xlon,vsize,'xlon','degrees_east', &
        'Longitude on Cross Points','longitude',lgetspace=.false.)
      call setup_var(var,xlat,vsize,'xlat','degrees_north', &
        'Latitude on Cross Points','latitude',lgetspace=.false.)
      call setup_var(var,mask,vsize,'mask','1', &
        'Land Mask','land_binary_mask',lgetspace=.false.)
      call setup_var(var,topo,vsize,'topo','m', &
        'Surface Model Elevation','surface_altitude',lgetspace=.false.)
      call setup_var(var,ps,vsize,'ps','hPa', &
        'Surface Pressure','surface_air_pressure',.true.,lgetspace=.false.)
      var(xlon)%rval => xlon_out
      var(xlat)%rval => xlat_out
      var(mask)%rval => mask_out
      var(topo)%rval => topo_out
      var(ps)%rval => ps_out
    else
      call setup_var(var,xlon,vsize,'xlon','degrees_east', &
        'Longitude on Cross Points','longitude')
      call setup_var(var,xlat,vsize,'xlat','degrees_north', &
        'Latitude on Cross Points','latitude')
      call setup_var(var,mask,vsize,'mask','1', &
        'Land Mask','land_binary_mask')
      call setup_var(var,topo,vsize,'topo','m', &
        'Surface Model Elevation','surface_altitude')
      call setup_var(var,ps,vsize,'ps','hPa', &
        'Surface Pressure','surface_air_pressure',.true.)
      xlon_out => var(xlon)%rval
      xlat_out => var(xlat)%rval
      mask_out => var(mask)%rval
      topo_out => var(topo)%rval
      ps_out => var(ps)%rval
    end if
  end subroutine setup_common_vars

  subroutine setup_var_2d(var,ivar,vsize,vname,vunit,long_name,standard_name, &
                         l_rec,cell_method,l_fill,rmissval,lgetspace)
    implicit none
    type(ncvariable2d_real) , dimension(:) , intent(inout) :: var
    integer , intent(in) :: ivar
    type(varspan) , intent(in) :: vsize
    character(len=*) , intent(in) :: vname , vunit , long_name , standard_name
    character(len=*) , intent(in) , optional :: cell_method
    logical , intent(in) , optional :: l_rec , l_fill , lgetspace
    real(rk4) , intent(in) , optional :: rmissval
    var(ivar)%vname = vname
    var(ivar)%vunit = vunit
    var(ivar)%long_name = long_name
    var(ivar)%standard_name = standard_name
    if ( present(cell_method) ) then
      var(ivar)%cell_method = cell_method
    end if
    if ( present(rmissval) .or. present(l_fill) ) then
      var(ivar)%lfillvalue = .true.
      if ( present(rmissval) ) then
        var(ivar)%rmissval = rmissval
      end if
    end if
    if ( present(l_rec) ) then
      var(ivar)%lrecords = l_rec
    end if
    ! This was the error on the IBM BlueGeneQ !
    ! The following WILL NOT WORK with IBM or PGI compiler.
    !
    ! if ( .not. present(a) .or (present(a) .and. a == .T.) )
    !
    if ( present(lgetspace) ) then
      if ( lgetspace ) then
        call getmem2d(var(ivar)%rval,vsize%j1,vsize%j2, &
              vsize%i1,vsize%i2,'ncout:setup_var:'//trim(var(ivar)%vname))
      end if
    else
      call getmem2d(var(ivar)%rval,vsize%j1,vsize%j2, &
            vsize%i1,vsize%i2,'ncout:setup_var:'//trim(var(ivar)%vname))
    end if
    var(ivar)%j1 = vsize%j1
    var(ivar)%j2 = vsize%j2
    var(ivar)%i1 = vsize%i1
    var(ivar)%i2 = vsize%i2
  end subroutine setup_var_2d

  subroutine setup_var_3d(var,ivar,vsize,vname,vunit,long_name,standard_name, &
                          l_rec,cell_method,l_fill,rmissval,lgetspace)
    implicit none
    type(ncvariable3d_real) , dimension(:) , intent(inout) :: var
    integer(ik4) , intent(in) :: ivar
    type(varspan) , intent(in) :: vsize
    character(len=*) , intent(in) :: vname , vunit , long_name , standard_name
    character(len=*) , intent(in) , optional :: cell_method
    logical , intent(in) , optional :: l_rec , l_fill , lgetspace
    real(rk4) , intent(in) , optional :: rmissval
    var(ivar)%vname = vname
    var(ivar)%vunit = vunit
    var(ivar)%long_name = long_name
    var(ivar)%standard_name = standard_name
    if ( present(cell_method) ) then
      var(ivar)%cell_method = cell_method
    end if
    if ( present(rmissval) .or. present(l_fill) ) then
      var(ivar)%lfillvalue = .true.
      if ( present(rmissval) ) var(ivar)%rmissval = rmissval
    end if
    if ( present(l_rec) ) then
      var(ivar)%lrecords = l_rec
    end if
    if ( present(lgetspace) ) then
      if ( lgetspace ) then
        call getmem3d(var(ivar)%rval,vsize%j1,vsize%j2,vsize%i1,vsize%i2, &
              vsize%k1,vsize%k2,'ncout:setup_var:'//trim(var(ivar)%vname))
      end if
    else
      call getmem3d(var(ivar)%rval,vsize%j1,vsize%j2,vsize%i1,vsize%i2, &
            vsize%k1,vsize%k2,'ncout:setup_var:'//trim(var(ivar)%vname))
    end if
    var(ivar)%j1 = vsize%j1
    var(ivar)%j2 = vsize%j2
    var(ivar)%i1 = vsize%i1
    var(ivar)%i2 = vsize%i2
    var(ivar)%k1 = vsize%k1
    var(ivar)%k2 = vsize%k2
  end subroutine setup_var_3d

  subroutine setup_var_4d(var,ivar,vsize,vname,vunit,long_name,standard_name, &
                          l_rec,cell_method,l_fill,rmissval,lgetspace)
    implicit none
    type(ncvariable4d_real) , dimension(:) , intent(inout) :: var
    integer(ik4) , intent(in) :: ivar
    type(varspan) , intent(in) :: vsize
    character(len=*) , intent(in) :: vname , vunit , long_name , standard_name
    character(len=*) , intent(in) , optional :: cell_method
    logical , intent(in) , optional :: l_rec , l_fill , lgetspace
    real(rk4) , intent(in) , optional :: rmissval
    var(ivar)%vname = vname
    var(ivar)%vunit = vunit
    var(ivar)%long_name = long_name
    var(ivar)%standard_name = standard_name
    if ( present(cell_method) ) then
      var(ivar)%cell_method = cell_method
    end if
    if ( present(rmissval) .or. present(l_fill) ) then
      var(ivar)%lfillvalue = .true.
      if ( present(rmissval) ) var(ivar)%rmissval = rmissval
    end if
    if ( present(l_rec) ) then
      var(ivar)%lrecords = l_rec
    end if
    if ( present(lgetspace) ) then
      if ( lgetspace ) then
        call getmem4d(var(ivar)%rval,vsize%j1,vsize%j2,vsize%i1,vsize%i2, &
              vsize%k1,vsize%k2,vsize%n1,vsize%n2,'ncout:setup_var:'//trim(var(ivar)%vname))
      end if
    else
      call getmem4d(var(ivar)%rval,vsize%j1,vsize%j2,vsize%i1,vsize%i2, &
            vsize%k1,vsize%k2,vsize%n1,vsize%n2,'ncout:setup_var:'//trim(var(ivar)%vname))
    end if
    var(ivar)%j1 = vsize%j1
    var(ivar)%j2 = vsize%j2
    var(ivar)%i1 = vsize%i1
    var(ivar)%i2 = vsize%i2
    var(ivar)%k1 = vsize%k1
    var(ivar)%k2 = vsize%k2
    var(ivar)%n1 = vsize%n1
    var(ivar)%n2 = vsize%n2
  end subroutine setup_var_4d

 subroutine dispose_output_streams
    implicit none
    integer(ik4) :: nstream , nfile
    if ( associated(v2dvar_atm) ) deallocate(v2dvar_atm)
    if ( associated(v3dvar_atm) ) deallocate(v3dvar_atm)
    if ( associated(v2dvar_srf) ) deallocate(v2dvar_srf)
    if ( associated(v3dvar_srf) ) deallocate(v3dvar_srf)
    if ( associated(v2dvar_sts) ) deallocate(v2dvar_sts)
    if ( associated(v3dvar_sts) ) deallocate(v3dvar_sts)
    if ( associated(v2dvar_rad) ) deallocate(v2dvar_rad)
    if ( associated(v3dvar_rad) ) deallocate(v3dvar_rad)
    if ( associated(v4dvar_rad) ) deallocate(v4dvar_rad)
    if ( associated(v2dvar_sub) ) deallocate(v2dvar_sub)
    if ( associated(v3dvar_sub) ) deallocate(v3dvar_sub)
    if ( associated(v2dvar_lak) ) deallocate(v2dvar_lak)
    if ( associated(v3dvar_lak) ) deallocate(v3dvar_lak)
    if ( associated(v2dvar_opt) ) deallocate(v2dvar_opt)
    if ( associated(v3dvar_opt) ) deallocate(v3dvar_opt)
    if ( associated(v2dvar_che) ) deallocate(v2dvar_che)
    if ( associated(v3dvar_che) ) deallocate(v3dvar_che)
    if ( associated(v2dvar_slaboc) ) deallocate(v2dvar_slaboc)
    do nstream = 1 , maxstreams
      do nfile = 1 , outstream(nstream)%nfiles
        call outstream_dispose(outstream(nstream)%ncout(nfile))
      end do
      deallocate(outstream(nstream)%ncout)
      deallocate(outstream(nstream)%cname_base)
      if ( associated(outstream(nstream)%ncvars%vlist) ) then
        deallocate(outstream(nstream)%ncvars%vlist)
      end if
    end do
    deallocate(outstream)
  end subroutine dispose_output_streams

  subroutine write_record_output_stream(istream,idate,ifile)
    implicit none
    integer(ik4) , intent(in) :: istream
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) , intent(in) , optional :: ifile
    real(rk8) , pointer , dimension(:,:) :: tmp2d
    real(rk8) , pointer , dimension(:,:,:) :: tmp3d
    real(rk8) , pointer , dimension(:,:,:,:) :: tmp4d
    real(rk8) , pointer , dimension(:,:) :: pnt2d => null( )
    real(rk8) , pointer , dimension(:,:,:) :: pnt3d => null( )
    real(rk8) , pointer , dimension(:,:,:,:) :: pnt4d => null( )
    class(ncvariable_standard) , pointer :: vp
    integer(ik4) :: ivar , jfile

    if ( .not. parallel_out .and. myid /= iocpu ) then
      do ivar = 1 , outstream(istream)%nvar
        vp => outstream(istream)%ncvars%vlist(ivar)%vp
        select type(vp)
          type is (ncvariable2d_real)
            if ( .not. vp%lrecords ) cycle
            call grid_collect(vp%rval,pnt2d, &
                              vp%j1,vp%j2,vp%i1,vp%i2,outstream(istream)%l_sub)
          type is (ncvariable3d_real)
            if ( .not. vp%lrecords ) cycle
            call grid_collect(vp%rval,pnt3d,vp%j1,vp%j2, &
                              vp%i1,vp%i2,vp%k1,vp%k2,outstream(istream)%l_sub)
          type is (ncvariable4d_real)
            if ( .not. vp%lrecords ) cycle
            call grid_collect(vp%rval,pnt4d,vp%j1,vp%j2, &
                              vp%i1,vp%i2,vp%k1,vp%k2,vp%n1,vp%n2)
          class default
            cycle
        end select
      end do

      ! If not parallel output, only the master proc writes output files
      return

    end if

    jfile = outstream(istream)%nfiles
    if ( present(ifile) ) jfile = ifile
    if ( jfile < 1 .or. jfile > outstream(istream)%nfiles ) then
      write (stderr,*) 'No such file in stream ',istream,' : ', ifile
      return
    end if

    call outstream_addrec(outstream(istream)%ncout(jfile),idate)

    if ( .not. parallel_out ) then
      if ( outstream(istream)%l_sub ) then
        pnt2d => io2dsg
        pnt3d => io3dsg
      else
        pnt2d => io2d
        pnt3d => io3d
        pnt4d => io4d
      end if
    end if

    do ivar = 1 , outstream(istream)%nvar

      vp => outstream(istream)%ncvars%vlist(ivar)%vp

      ! If not parallel output, collect data

      if ( .not. parallel_out ) then
        select type(vp)
          type is (ncvariable2d_real)
            if ( .not. vp%lrecords ) cycle
            call grid_collect(vp%rval,pnt2d, &
                              vp%j1,vp%j2,vp%i1,vp%i2,outstream(istream)%l_sub)
            vp%j1 = outstream(istream)%jg1
            vp%j2 = outstream(istream)%jg2
            vp%i1 = outstream(istream)%ig1
            vp%i2 = outstream(istream)%ig2
            tmp2d => vp%rval
            vp%rval => pnt2d
          type is (ncvariable3d_real)
            if ( .not. vp%lrecords ) cycle
            call grid_collect(vp%rval,pnt3d,vp%j1,vp%j2, &
                              vp%i1,vp%i2,vp%k1,vp%k2,outstream(istream)%l_sub)
            vp%j1 = outstream(istream)%jg1
            vp%j2 = outstream(istream)%jg2
            vp%i1 = outstream(istream)%ig1
            vp%i2 = outstream(istream)%ig2
            tmp3d => vp%rval
            vp%rval => pnt3d
          type is (ncvariable4d_real)
            if ( .not. vp%lrecords ) cycle
            call grid_collect(vp%rval,pnt4d,vp%j1,vp%j2, &
                              vp%i1,vp%i2,vp%k1,vp%k2,vp%n1,vp%n2)
            vp%j1 = outstream(istream)%jg1
            vp%j2 = outstream(istream)%jg2
            vp%i1 = outstream(istream)%ig1
            vp%i2 = outstream(istream)%ig2
            tmp4d => vp%rval
            vp%rval => pnt4d
          class default
            cycle
        end select
      else
        select type(vp)
          type is (ncvariable2d_real)
            if ( .not. vp%lrecords ) cycle
          type is (ncvariable3d_real)
            if ( .not. vp%lrecords ) cycle
          type is (ncvariable4d_real)
            if ( .not. vp%lrecords ) cycle
          class default
            cycle
        end select
      end if

#ifdef DEBUG
    if ( debug_level > 2 ) then
      write(ndebug+myid,*) &
              'Writing var ',trim(vp%vname),' at time ',tochar(idate)
      select type(vp)
        type is (ncvariable2d_real)
          write(ndebug+myid,*) 'Max Value : ', maxval(vp%rval)
          write(ndebug+myid,*) 'Min Value : ', minval(vp%rval)
        type is (ncvariable3d_real)
          write(ndebug+myid,*) 'Max Value : ', maxval(vp%rval)
          write(ndebug+myid,*) 'Min Value : ', minval(vp%rval)
        type is (ncvariable4d_real)
          write(ndebug+myid,*) 'Max Value : ', maxval(vp%rval)
          write(ndebug+myid,*) 'Min Value : ', minval(vp%rval)
        class default
      end select
    end if
#endif
      call outstream_writevar(outstream(istream)%ncout(jfile),vp)

      ! Reset pointers

      if ( .not. parallel_out ) then
        select type(vp)
          type is (ncvariable2d_real)
            vp%rval => tmp2d
            vp%j1 = outstream(istream)%jl1
            vp%j2 = outstream(istream)%jl2
            vp%i1 = outstream(istream)%il1
            vp%i2 = outstream(istream)%il2
          type is (ncvariable3d_real)
            vp%rval => tmp3d
            vp%j1 = outstream(istream)%jl1
            vp%j2 = outstream(istream)%jl2
            vp%i1 = outstream(istream)%il1
            vp%i2 = outstream(istream)%il2
          type is (ncvariable4d_real)
            vp%rval => tmp4d
            vp%j1 = outstream(istream)%jl1
            vp%j2 = outstream(istream)%jl2
            vp%i1 = outstream(istream)%il1
            vp%i2 = outstream(istream)%il2
          class default
            cycle
        end select
      end if

    end do

  end subroutine write_record_output_stream

  subroutine writevar2d_output_stream(istream,vp,ifile)
    implicit none
    integer(ik4) , intent(in) :: istream
    integer(ik4) , intent(in) , optional :: ifile
    type(ncvariable2d_real) , intent(inout) :: vp
    real(rk8) , pointer , dimension(:,:) :: tmp2d
    real(rk8) , pointer , dimension(:,:) :: pnt2d => null( )
    integer(ik4) :: jfile

    if ( .not. parallel_out .and. myid /= iocpu ) then
      call grid_collect(vp%rval,pnt2d, &
                        vp%j1,vp%j2,vp%i1,vp%i2,outstream(istream)%l_sub)
      ! If not parallel output, only the master proc writes output files
      return
    end if

    jfile = outstream(istream)%nfiles
    if ( present(ifile) ) jfile = ifile
    if ( jfile < 1 .or. jfile > outstream(istream)%nfiles ) then
      write (stderr,*) 'No such file in stream ',istream,' : ', ifile
      return
    end if

    if ( .not. parallel_out ) then
      if ( outstream(istream)%l_sub ) then
        pnt2d => io2dsg
      else
        pnt2d => io2d
      end if
    end if

    ! If not parallel output, collect data

    if ( .not. parallel_out ) then
      call grid_collect(vp%rval,pnt2d, &
                        vp%j1,vp%j2,vp%i1,vp%i2,outstream(istream)%l_sub)
      vp%j1 = outstream(istream)%jg1
      vp%j2 = outstream(istream)%jg2
      vp%i1 = outstream(istream)%ig1
      vp%i2 = outstream(istream)%ig2
      tmp2d => vp%rval
      vp%rval => pnt2d
    end if

#ifdef DEBUG
    if ( debug_level > 2 ) then
      write(ndebug+myid,*) 'Writing var ',trim(vp%vname)
    end if
#endif
    call outstream_writevar(outstream(istream)%ncout(jfile),vp)

    ! Reset pointers

    if ( .not. parallel_out ) then
      vp%rval => tmp2d
      vp%j1 = outstream(istream)%jl1
      vp%j2 = outstream(istream)%jl2
      vp%i1 = outstream(istream)%il1
      vp%i2 = outstream(istream)%il2
    end if

  end subroutine writevar2d_output_stream

  subroutine writevar3d_output_stream(istream,vp,ifile)
    implicit none
    integer(ik4) , intent(in) :: istream
    integer(ik4) , intent(in) , optional :: ifile
    type(ncvariable3d_real) , intent(inout) :: vp
    real(rk8) , pointer , dimension(:,:,:) :: tmp3d
    real(rk8) , pointer , dimension(:,:,:) :: pnt3d => null( )
    integer(ik4) :: jfile

    if ( .not. parallel_out .and. myid /= iocpu ) then
      call grid_collect(vp%rval,pnt3d,vp%j1,vp%j2, &
                        vp%i1,vp%i2,vp%k1,vp%k2,outstream(istream)%l_sub)
      ! If not parallel output, only the master proc writes output files
      return
    end if

    jfile = outstream(istream)%nfiles
    if ( present(ifile) ) jfile = ifile
    if ( jfile < 1 .or. jfile > outstream(istream)%nfiles ) then
      write (stderr,*) 'No such file in stream ',istream,' : ', ifile
      return
    end if

    if ( .not. parallel_out ) then
      if ( outstream(istream)%l_sub ) then
        pnt3d => io3dsg
      else
        pnt3d => io3d
      end if
    end if

    ! If not parallel output, collect data

    if ( .not. parallel_out ) then
      call grid_collect(vp%rval,pnt3d,vp%j1,vp%j2, &
                        vp%i1,vp%i2,vp%k1,vp%k2,outstream(istream)%l_sub)
      vp%j1 = outstream(istream)%jg1
      vp%j2 = outstream(istream)%jg2
      vp%i1 = outstream(istream)%ig1
      vp%i2 = outstream(istream)%ig2
      tmp3d => vp%rval
      vp%rval => pnt3d
    end if

#ifdef DEBUG
    if ( debug_level > 2 ) then
      write(ndebug+myid,*) 'Writing var ',trim(vp%vname)
    end if
#endif
    call outstream_writevar(outstream(istream)%ncout(jfile),vp)

    ! Reset pointers

    if ( .not. parallel_out ) then
      vp%rval => tmp3d
      vp%j1 = outstream(istream)%jl1
      vp%j2 = outstream(istream)%jl2
      vp%i1 = outstream(istream)%il1
      vp%i2 = outstream(istream)%il2
    end if
  end subroutine writevar3d_output_stream

  subroutine writevar4d_output_stream(istream,vp,ifile)
    implicit none
    integer(ik4) , intent(in) :: istream
    integer(ik4) , intent(in) , optional :: ifile
    type(ncvariable4d_real) , intent(inout) :: vp
    real(rk8) , pointer , dimension(:,:,:,:) :: tmp4d
    real(rk8) , pointer , dimension(:,:,:,:) :: pnt4d => null( )
    integer(ik4) :: jfile

    if ( .not. parallel_out .and. myid /= iocpu ) then
      call grid_collect(vp%rval,pnt4d,vp%j1,vp%j2, &
                        vp%i1,vp%i2,vp%k1,vp%k2,vp%n1,vp%n2)
      ! If not parallel output, only the master proc writes output files
      return
    end if

    jfile = outstream(istream)%nfiles
    if ( present(ifile) ) jfile = ifile
    if ( jfile < 1 .or. jfile > outstream(istream)%nfiles ) then
      write (stderr,*) 'No such file in stream ',istream,' : ', ifile
      return
    end if

    if ( .not. parallel_out ) then
      pnt4d => io4d
    end if

    ! If not parallel output, collect data

    if ( .not. parallel_out ) then
      call grid_collect(vp%rval,pnt4d,vp%j1,vp%j2, &
                        vp%i1,vp%i2,vp%k1,vp%k2,vp%n1,vp%n2)
      vp%j1 = outstream(istream)%jg1
      vp%j2 = outstream(istream)%jg2
      vp%i1 = outstream(istream)%ig1
      vp%i2 = outstream(istream)%ig2
      tmp4d => vp%rval
      vp%rval => pnt4d
    end if

#ifdef DEBUG
    if ( debug_level > 2 ) then
      write(ndebug+myid,*) 'Writing var ',trim(vp%vname)
    end if
#endif
    call outstream_writevar(outstream(istream)%ncout(jfile),vp)

    ! Reset pointers

    if ( .not. parallel_out ) then
      vp%rval => tmp4d
      vp%j1 = outstream(istream)%jl1
      vp%j2 = outstream(istream)%jl2
      vp%i1 = outstream(istream)%il1
      vp%i2 = outstream(istream)%il2
    end if
  end subroutine writevar4d_output_stream

end module mod_ncout
