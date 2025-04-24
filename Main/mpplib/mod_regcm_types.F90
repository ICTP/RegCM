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

module mod_regcm_types
  use mod_realkinds
  use mod_intkinds

  implicit none

  public

!
! Storage for all the 3d prognostic variables in two
! timesteps and all the 2d variables and constants
!
  type masked_comm
    integer(ik4) :: linear_communicator
    logical, pointer, contiguous, dimension(:,:) :: gmask => null( )
    logical, pointer, contiguous, dimension(:,:,:) :: sgmask => null( )
    logical, pointer, contiguous, dimension(:,:) :: global_gmask => null( )
    logical, pointer, contiguous, dimension(:,:,:) :: global_sgmask => null( )
    logical, pointer, contiguous, dimension(:,:) :: global_out_sgmask => null( )
    integer(ik4), pointer, contiguous, dimension(:) :: linear_npoint_g => null( )
    integer(ik4), pointer, contiguous, dimension(:) :: linear_displ_g => null( )
    integer(ik4), pointer, contiguous, dimension(:) :: cartesian_npoint_g => null( )
    integer(ik4), pointer, contiguous, dimension(:) :: cartesian_displ_g => null( )
    integer(ik4), pointer, contiguous, dimension(:) :: linear_npoint_sg => null( )
    integer(ik4), pointer, contiguous, dimension(:) :: linear_displ_sg => null( )
    integer(ik4), pointer, contiguous, dimension(:) :: cartesian_npoint_sg => null( )
    integer(ik4), pointer, contiguous, dimension(:) :: cartesian_displ_sg => null( )
  end type masked_comm

  type model_area
    logical :: bandflag
    logical :: crmflag
    logical :: has_bdy
    logical :: has_bdyleft, has_bdyright, has_bdytop, has_bdybottom
    logical :: has_bdytopleft, has_bdytopright
    logical :: has_bdybottomleft, has_bdybottomright
    integer(ik4), dimension(2) :: location
    integer(ik4) :: left, right, top, bottom
    integer(ik4) :: topleft, topright, bottomleft, bottomright
    integer(ik4) :: ibt1, ibt2, ibt3, ibt4, ibt6
    integer(ik4) :: ibb1, ibb2, ibb3, ibb4, ibb6
    integer(ik4) :: jbl1, jbl2, jbl3, jbl4, jbl6
    integer(ik4) :: jbr1, jbr2, jbr3, jbr4, jbr6
  end type model_area

  type domain
    real(rkx), pointer, contiguous, dimension(:,:) :: ht => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: lndcat => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: lndtex => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: mask => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: area => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: dlat => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: dlon => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ulat => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ulon => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: vlat => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: vlon => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: xlat => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: xlon => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: msfu => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: msfv => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: hx => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: hy => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: msfx => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: msfd => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: coriol => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: coriou => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: coriov => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ef => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ddx => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ddy => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ex => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: crx => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: cry => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: dmdy => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: dmdx => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: xmsf => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: dmsf => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: snowam => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: smoist => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rmoist => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rts => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: dhlake => null( )
    integer(ik4), pointer, contiguous, dimension(:,:) :: ldmsk => null( )
    integer(ik4), pointer, contiguous, dimension(:,:) :: iveg => null( )
    integer(ik4), pointer, contiguous, dimension(:,:) :: itex => null( )
  end type domain

  type domain_subgrid
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ht => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: lndcat => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: lndtex => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: xlat => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: xlon => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: mask => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: area => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dhlake => null( )
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: ldmsk => null( )
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: iveg => null( )
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: itex => null( )
  end type domain_subgrid

  type mass_divergence
    real(rkx), pointer, contiguous, dimension(:,:,:) :: cr => null( ) ! cross points
  end type mass_divergence

  type atmosphere
    real(rkx), pointer, contiguous, dimension(:,:,:) :: u => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: v => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ux => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vx => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: w => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pai => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: t => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: p => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rho => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tvirt => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tetav => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tke => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qs => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qx => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: trac => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: zeta => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: zetaf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dz => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: fmz => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: fmzf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tten => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: uten => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vten => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tketen => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qxten => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chiten => null( )
  end type atmosphere

  type reference_atmosphere
    real(rkx), pointer, contiguous, dimension(:,:) :: ps => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: psdot => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: t => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pr => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rho => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: z => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: zd => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rhof => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: zf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dzf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dprddx => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dprddy => null( )
  end type reference_atmosphere

  type atmstate_a
    real(rkx), pointer, contiguous, dimension(:,:,:) :: u => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: v => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: w => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: t => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tke => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pp => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pr => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rho => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qx => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chi => null( )
  end type atmstate_a

  type atmstate_b
    real(rkx), pointer, contiguous, dimension(:,:,:) :: u => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: v => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: w => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: t => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tke => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pp => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pr => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qx => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chi => null( )
  end type atmstate_b

  type atmstate_c
    real(rkx), pointer, contiguous, dimension(:,:,:) :: u => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: v => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: w => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: t => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tke => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pp => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qx => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chi => null( )
  end type atmstate_c

  type atmstate_tendency
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: u => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: v => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: w => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: t => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: tke => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: pp => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:,:) :: qx => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:,:) :: chi => null( )
  end type atmstate_tendency

  type crosswind_tendency
    real(rkx), pointer, contiguous, dimension(:,:,:) :: u => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: v => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ud => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vd => null( )
  end type crosswind_tendency

  type atmstate_decoupled
    real(rkx), pointer, contiguous, dimension(:,:,:) :: uc => null( )  ! pressure
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vc => null( )  ! pressure
    real(rkx), pointer, contiguous, dimension(:,:,:) :: umc => null( ) ! mapfactor
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vmc => null( ) ! mapfactor
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ud => null( )  ! De-coupled
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vd => null( )  ! De-coupled
    real(rkx), pointer, contiguous, dimension(:,:,:) :: umd => null( ) ! mapfactor
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vmd => null( ) ! mapfactor
    real(rkx), pointer, contiguous, dimension(:,:,:) :: w => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: t => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tv => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qx => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pp => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pr => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rho => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chi => null( )
  end type atmstate_decoupled

  type tcm_state
    ! Vertical momentum diffusivity
    real(rkx), pointer, contiguous, dimension(:,:,:) :: kzm => null( )    ! (m^2/s)
    ! Vertical scalar diffusivity
    real(rkx), pointer, contiguous, dimension(:,:,:) :: kth => null( )    ! (m^2/s)
  end type tcm_state

  type tendiag
    real(rkx), pointer, contiguous, dimension(:,:,:) :: adh => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: adv => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tbl => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dif => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: bdy => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: con => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: adi => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rad => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: lsc => null( )
  end type tendiag

  type qendiag
    real(rkx), pointer, contiguous, dimension(:,:,:) :: adh => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: adv => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tbl => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dif => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: bdy => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: con => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: adi => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rad => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: lsc => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qcl => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qcr => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: acr => null( )
  end type qendiag

  type surfstate
    real(rkx), pointer, contiguous, dimension(:,:) :: psa => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: psb => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: psc => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: psdota => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: psdotb => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: tg => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: rainc => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: rainnc => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: snownc => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: grplnc => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: hailnc => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: hfx => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: qfx => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: tgbb => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: q2m => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: uvdrag => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ustar => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: u10m => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: v10m => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: w10m => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: zo => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ram1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: rah1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: br => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: uz0 => null( )    ! MYJ SF layer
    real(rkx), pointer, contiguous, dimension(:,:) :: vz0 => null( )    ! MYJ SF layer
    real(rkx), pointer, contiguous, dimension(:,:) :: thz0 => null( )   ! MYJ SF layer
    real(rkx), pointer, contiguous, dimension(:,:) :: qz0 => null( )    ! MYJ SF layer
    real(rkx), pointer, contiguous, dimension(:,:) :: dtrnof => null( ) ! rnoff
  end type surfstate

  type slice
    real(rkx), pointer, contiguous, dimension(:,:,:) :: th3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: th700 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tb3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tv3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pb3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pf3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rhob3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ubx3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vbx3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: wb3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ppb3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: wpx3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ubd3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vbd3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rhb3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qsb3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qxb3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: zq => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: za => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dzq => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: rhox2d => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: tp2d => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ps2d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chib3d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tkepbl => null( )
  end type slice

  type v3dbound
    real(rkx), pointer, contiguous, dimension(:,:,:) :: b0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: b1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: bt => null( )
  end type v3dbound

  type v2dbound
    real(rkx), pointer, contiguous, dimension(:,:) :: b0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: b1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: bt => null( )
  end type v2dbound

  type nhboundhelp
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tvirt => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ps => null( )
  end type nhboundhelp

  type bound_area
    logical :: havebound
    logical, pointer, contiguous, dimension(:,:) :: bsouth => null( )
    logical, pointer, contiguous, dimension(:,:) :: bnorth => null( )
    logical, pointer, contiguous, dimension(:,:) :: beast => null( )
    logical, pointer, contiguous, dimension(:,:) :: bwest => null( )
    integer(ik4) :: ns, nn, ne, nw
    integer(ik4) :: nsp
    integer(ik4), pointer, contiguous, dimension(:,:) :: ibnd => null( )
  end type bound_area

  type exp_data
    real(rk8), pointer, contiguous, dimension(:,:) :: psfc => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: tsfc => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: qsfc => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: swrd => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: lwrd => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: dlwr => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: lhfx => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: shfx => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: prec => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: wndu => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: wndv => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: rnof => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: snof => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: taux => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: tauy => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: wspd => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: wdir => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: ustr => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: nflx => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: sflx => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: snow => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: dswr => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: rhoa => null( )
  end type exp_data

  type exp_data3d
    real(rk8), pointer, contiguous, dimension(:,:,:) :: t => null( )
    real(rk8), pointer, contiguous, dimension(:,:,:) :: q => null( )
    real(rk8), pointer, contiguous, dimension(:,:,:) :: u => null( )
    real(rk8), pointer, contiguous, dimension(:,:,:) :: v => null( )
    real(rk8), pointer, contiguous, dimension(:,:,:) :: w => null( )
    real(rk8), pointer, contiguous, dimension(:,:,:) :: cldfrc => null( )
    real(rk8), pointer, contiguous, dimension(:,:,:) :: cldlwc => null( )
  end type exp_data3d

  type imp_data
    real(rkx), pointer, contiguous, dimension(:,:) :: sst => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: sit => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: msk => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: zo => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ustar => null( )
  end type imp_data

  type lm_state
    real(rkx), pointer, contiguous, dimension(:,:,:) :: gwet => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: sw => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ssw => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rsw => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tsw => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: lncl => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ldew => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tgbb => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tgrd => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tgbrd => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: taf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tlef => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: sfice => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: snag => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: sncv => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: scvk => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: emisv => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: sent => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: evpr => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: deltat => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: deltaq => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: drag => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: prcp => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: snwm => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: trnof => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: sigf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: sfcp => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: srnof => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: xlai => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: q2m => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: t2m => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: u10m => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: v10m => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ram1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rah1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: br => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: taux => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tauy => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: swalb => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: lwalb => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: urlwf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: swdiralb => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: lwdiralb => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: swdifalb => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: lwdifalb => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: wt => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: eta => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: hi => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: hsnow => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: um10 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: tlake => null( )
    logical, pointer, contiguous, dimension(:,:,:) :: lakmsk => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: deltas => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tdeltas => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tskin => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: sst => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: zo => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ustar => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: w10m => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rhoa => null( )
#ifdef CLM45
    real(rkx), pointer, contiguous, dimension(:,:,:) :: hfso => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: vocemiss => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: dustemiss => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: ddepv => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: sw_vol => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: tsoi => null( )
#endif
  end type lm_state

  type lm_exchange
    real(rkx), pointer, contiguous, dimension(:,:) :: ssw2da => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: sfracv2d => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: sfracb2d => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: sfracs2d => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: svegfrac2d => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: sxlai2d => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: sw_vol => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tsoi => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: xlat => null( )     ! mddom%xlat
    real(rkx), pointer, contiguous, dimension(:,:) :: xlon => null( )     ! mddom%xlon
    real(rkx), pointer, contiguous, dimension(:,:) :: lndcat => null( )   ! mddom%lndcat
    real(rkx), pointer, contiguous, dimension(:,:) :: ht => null( )       ! mddom%ht
    real(rkx), pointer, contiguous, dimension(:,:) :: snowam => null( )   ! mddom%snowam
    real(rkx), pointer, contiguous, dimension(:,:) :: smoist => null( )   ! mddom%smoist
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rmoist => null( ) ! mddom%rmoist
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rts => null( )    ! mddom%rts
    integer(ik4), pointer, contiguous, dimension(:,:) :: iveg => null( )  ! mddom%iveg
    integer(ik4), pointer, contiguous, dimension(:,:) :: itex => null( )  ! mddom%itex
    integer(ik4), pointer, contiguous, dimension(:,:) :: ldmsk => null( ) ! mddom%ldmsk
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ht1 => null( )    ! mdsub%ht
    real(rkx), pointer, contiguous, dimension(:,:,:) :: lndcat1 => null( ) ! mdsub%lndcat
    real(rkx), pointer, contiguous, dimension(:,:,:) :: area1 => null( )   ! mdsub%area
    real(rkx), pointer, contiguous, dimension(:,:,:) :: xlat1 => null( )   ! mdsub%xlat
    real(rkx), pointer, contiguous, dimension(:,:,:) :: xlon1 => null( )   ! mdsub%xlon
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dhlake1 => null( ) ! mdsub%dhlake
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: ldmsk1 => null( ) ! mdsub%ldmsk
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: iveg1 => null( )  ! mdsub%iveg
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: itex1 => null( )  ! mdsub%itex
    integer(ik4), pointer, contiguous, dimension(:,:) :: icplmsk => null( )  ! cplmsk
    real(rkx), pointer, contiguous, dimension(:,:) :: patm => null( ) ! atms%pb3d(kz)
    real(rkx), pointer, contiguous, dimension(:,:) :: uatm => null( ) ! atms%ubx3d(kz)
    real(rkx), pointer, contiguous, dimension(:,:) :: vatm => null( ) ! atms%vbx3d(kz)
    real(rkx), pointer, contiguous, dimension(:,:) :: tatm => null( ) ! atms%tb3d(kz)
    real(rkx), pointer, contiguous, dimension(:,:) :: thatm => null( )! atms%th3d(kz)
    real(rkx), pointer, contiguous, dimension(:,:) :: qvatm => null( )! atms%qxb3d(kz,iqv)
    real(rkx), pointer, contiguous, dimension(:,:) :: hgt => null( )  ! za(kz)
    real(rkx), pointer, contiguous, dimension(:,:) :: hpbl => null( ) ! zpbl
    real(rkx), pointer, contiguous, dimension(:,:) :: hfx => null( )  ! sfs%hfx
    real(rkx), pointer, contiguous, dimension(:,:) :: qfx => null( )  ! sfs%qfx
    real(rkx), pointer, contiguous, dimension(:,:) :: totcf => null( ) ! totcf
    real(rkx), pointer, contiguous, dimension(:,:) :: sfps => null( )  ! sfs%psb
    real(rkx), pointer, contiguous, dimension(:,:) :: sfta => null( )  ! atms%ts2d
    real(rkx), pointer, contiguous, dimension(:,:) :: uvdrag => null( ) ! sfs%uvdrag
    real(rkx), pointer, contiguous, dimension(:,:) :: tg => null( )     ! sfs%tg
    real(rkx), pointer, contiguous, dimension(:,:) :: tgbb => null( )   ! sfs%tgbb
    real(rkx), pointer, contiguous, dimension(:,:) :: ram1 => null( )   ! sfs%ram1
    real(rkx), pointer, contiguous, dimension(:,:) :: rah1 => null( )   ! sfs%rah1
    real(rkx), pointer, contiguous, dimension(:,:) :: br => null( )     ! sfs%br
    real(rkx), pointer, contiguous, dimension(:,:) :: u10m => null( )   ! sfs%u10m
    real(rkx), pointer, contiguous, dimension(:,:) :: v10m => null( )   ! sfs%v10m
    real(rkx), pointer, contiguous, dimension(:,:) :: q2m => null( )    ! sfs%q2m
    real(rkx), pointer, contiguous, dimension(:,:) :: rhox => null( )   ! rhox2d
    real(rkx), pointer, contiguous, dimension(:,:) :: rswf => null( )   ! fsw
    real(rkx), pointer, contiguous, dimension(:,:) :: rlwf => null( )   ! flw
    real(rkx), pointer, contiguous, dimension(:,:) :: dwrlwf => null( ) ! flwd
    real(rkx), pointer, contiguous, dimension(:,:) :: zencos => null( )   ! coszrs
    real(rkx), pointer, contiguous, dimension(:,:) :: dtrnof => null( ) ! rnoff
    real(rkx), pointer, contiguous, dimension(:,:) :: ncprate => null( )  ! pptnc
    real(rkx), pointer, contiguous, dimension(:,:) :: cprate => null( )   ! cprate
    real(rkx), pointer, contiguous, dimension(:,:) :: snwrat => null( )   ! snwrat
    real(rkx), pointer, contiguous, dimension(:,:) :: csrate => null( )   ! csrate
    real(rkx), pointer, contiguous, dimension(:,:) :: grprat => null( )   ! grprat
    real(rkx), pointer, contiguous, dimension(:,:) :: hairat => null( )   ! hairat
    real(rkx), pointer, contiguous, dimension(:,:) :: vegswab => null( )  ! sabveg
    real(rkx), pointer, contiguous, dimension(:,:) :: lwalb => null( )    ! albvl
    real(rkx), pointer, contiguous, dimension(:,:) :: swalb => null( )    ! albvs
    real(rkx), pointer, contiguous, dimension(:,:) :: swdiralb => null( ) ! aldirs
    real(rkx), pointer, contiguous, dimension(:,:) :: swdifalb => null( ) ! aldifs
    real(rkx), pointer, contiguous, dimension(:,:) :: lwdiralb => null( ) ! aldirl
    real(rkx), pointer, contiguous, dimension(:,:) :: lwdifalb => null( ) ! aldifl
    real(rkx), pointer, contiguous, dimension(:,:) :: swdir => null( )    ! solvs
    real(rkx), pointer, contiguous, dimension(:,:) :: swdif => null( )    ! solvsd
    real(rkx), pointer, contiguous, dimension(:,:) :: lwdir => null( )    ! solvl
    real(rkx), pointer, contiguous, dimension(:,:) :: lwdif => null( )    ! solvld
    real(rkx), pointer, contiguous, dimension(:,:) :: solinc => null( )   ! sinc
    real(rkx), pointer, contiguous, dimension(:,:) :: solar => null( )    ! solis
    real(rkx), pointer, contiguous, dimension(:,:) :: emissivity => null( ) ! emiss
    real(rkx), pointer, contiguous, dimension(:,:) :: deltaq => null( )     ! sdelq
    real(rkx), pointer, contiguous, dimension(:,:) :: deltat => null( )     ! sdelt
    real(rkx), pointer, contiguous, dimension(:,:,:) :: drydepflx => null( ) ! drydepflx
    real(rkx), pointer, contiguous, dimension(:,:,:) :: wetdepflx => null( ) ! wetdepflx
    integer(ik4), pointer, contiguous, dimension(:) :: idust => null( ) ! dust indices
    real(rkx), pointer, contiguous, dimension(:,:) :: zo => null( )     ! zo
    real(rkx), pointer, contiguous, dimension(:,:) :: ustar => null( )  ! ustar
    real(rkx), pointer, contiguous, dimension(:,:) :: w10m => null( )   ! w10m
#ifdef CLM
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dep_vels => null( )
#endif
  end type lm_exchange

  type mod_2_rad
    real(rkx), pointer, contiguous, dimension(:,:,:) :: cldfrc => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: cldlwc => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tatms => null( )     ! atms%tb3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rhatms => null( )    ! atms%rhb3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rhoatms => null( )   ! atms%rhob3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: phatms => null( )    ! atms%pb3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pfatms => null( )    ! atms%pf3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: za => null( )        ! atms%za
    real(rkx), pointer, contiguous, dimension(:,:,:) :: zq => null( )        ! atms%zq
    real(rkx), pointer, contiguous, dimension(:,:,:) :: deltaz => null( )    ! atms%dzq
    real(rkx), pointer, contiguous, dimension(:,:) :: psatms => null( )      ! atms%ps2d
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qxatms => null( )  ! atms%qxb3d
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chiatms => null( ) ! atms%chib3d
    real(rkx), pointer, contiguous, dimension(:,:) :: tg => null( )     ! sfs%tgbb
    real(rkx), pointer, contiguous, dimension(:,:) :: xlat => null( )   ! mddom%xlat
    real(rkx), pointer, contiguous, dimension(:,:) :: xlon => null( )   ! mddom%xlon
    real(rkx), pointer, contiguous, dimension(:,:) :: ht => null( )     ! mddom%ht
    real(rkx), pointer, contiguous, dimension(:,:) :: ptrop => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: coszrs => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: albvs => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: albvl => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: aldirs => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: aldifs => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: aldirl => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: aldifl => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: emiss => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: ps0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: bps0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: btv0 => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: bps1 => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: btv1 => null( )
    integer(ik4), pointer, contiguous, dimension(:,:) :: ldmsk => null( )
  end type mod_2_rad

  type rad_2_mod
    real(rkx), pointer, contiguous, dimension(:,:) :: solis => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: sinc => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: sabveg => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: sols => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: soll => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: solvs => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: solvsd => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: solvl => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: solvld => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: fsw => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: flw => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: flwd => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: totcf => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: heatrt => null( )
  end type rad_2_mod

  type mod_2_cum
    real(rkx), pointer, contiguous, dimension(:,:) :: ht => null( )        ! mddom%ht
    real(rkx), pointer, contiguous, dimension(:,:) :: psa => null( )       ! sfs%psa
    real(rkx), pointer, contiguous, dimension(:,:) :: psb => null( )       ! sfs%psb
    real(rkx), pointer, contiguous, dimension(:,:) :: psdotb => null( )    ! sfs%psdotb
    real(rkx), pointer, contiguous, dimension(:,:) :: psf => null( )       ! atms%ps2d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pas => null( )     ! atms%pb3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pasf => null( )    ! atms%pf3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: zas => null( )     ! atms%za
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tas => null( )     ! atms%tb3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: uas => null( )     ! atms%ubx3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vas => null( )     ! atms%vbx3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: was => null( )     ! atms%wx3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: wpas => null( )    ! atms%wpx3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qsas => null( )    ! atms%qsb3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tkeas => null( )   ! atms%tke
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rhoas => null( )   ! atms%rhob3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: zfs => null( )     ! atms%zq
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dzq => null( )     ! atms%dzq
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qq1 => null( )     ! atm1%q
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qdot => null( )    ! qdot
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qxas => null( )  ! atms%qxb3d
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chias => null( ) ! atms%chib3d
    real(rkx), pointer, contiguous, dimension(:,:) :: qfx => null( )       ! sfs%qfx
    real(rkx), pointer, contiguous, dimension(:,:) :: hfx => null( )       ! sfs%hfx
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tten => null( )     ! aten%t
    real(rkx), pointer, contiguous, dimension(:,:,:) :: uten => null( )     ! aten%u
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vten => null( )     ! aten%v
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qxten => null( )  ! aten%qx
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: dynqx => null( )  ! aten%qx
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chiten => null( ) ! aten%chi
    real(rkx), pointer, contiguous, dimension(:,:,:) :: heatrt => null( )   ! radiation
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ccn => null( )      ! ccn
    integer(ik4), pointer, contiguous, dimension(:,:) :: ktrop => null( )
    integer(ik4), pointer, contiguous, dimension(:,:) :: ldmsk => null( )   ! mddom
  end type mod_2_cum

  type cum_2_mod
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tten => null( )     ! aten%t
    real(rkx), pointer, contiguous, dimension(:,:,:) :: uten => null( )     ! aten%u
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vten => null( )     ! aten%v
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qxten => null( )  ! aten%qx
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chiten => null( ) ! aten%chi
    real(rkx), pointer, contiguous, dimension(:,:) :: rainc => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: pcratec => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: sratec => null( )
    real(rkx), pointer, contiguous, dimension(:,:) :: trrate => null( )     ! trrate
    real(rkx), pointer, contiguous, dimension(:,:,:) :: convpr => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: cldfrc => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: cldlwc => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: q_detr => null( )
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rain_cc => null( )
    integer(ik4), pointer, contiguous, dimension(:,:) :: kcumtop => null( )
    integer(ik4), pointer, contiguous, dimension(:,:) :: kcumbot => null( )
  end type cum_2_mod

  type mod_2_pbl
    real(rkx), pointer, contiguous, dimension(:,:) :: coriol => null( )    ! mddom%coriol
    integer(ik4), pointer, contiguous, dimension(:,:) :: ldmsk => null( )  ! mddom%ldmsk
    real(rkx), pointer, contiguous, dimension(:,:) :: ht => null( )        ! mddom%ht
    real(rkx), pointer, contiguous, dimension(:,:) :: psb => null( )       ! sfs%psb
    real(rkx), pointer, contiguous, dimension(:,:) :: psdotb => null( )    ! sfs%psdotb
    real(rkx), pointer, contiguous, dimension(:,:) :: tg => null( )        ! sfs%tgbb
    real(rkx), pointer, contiguous, dimension(:,:) :: q2m => null( )       ! sfs%q2m
    real(rkx), pointer, contiguous, dimension(:,:) :: u10m => null( )      ! sfs%u10m
    real(rkx), pointer, contiguous, dimension(:,:) :: v10m => null( )      ! sfs%v10m
    real(rkx), pointer, contiguous, dimension(:,:) :: qfx => null( )       ! sfs%qfx
    real(rkx), pointer, contiguous, dimension(:,:) :: hfx => null( )       ! sfs%hfx
    real(rkx), pointer, contiguous, dimension(:,:) :: uvdrag => null( )    ! sfs%uvdrag
    real(rkx), pointer, contiguous, dimension(:,:) :: zo => null( )        ! sfs%zo
    real(rkx), pointer, contiguous, dimension(:,:) :: ustar => null( )     ! sfs%ustar
    real(rkx), pointer, contiguous, dimension(:,:) :: ram1 => null( )      ! sfs%ram1
    real(rkx), pointer, contiguous, dimension(:,:) :: rah1 => null( )      ! sfs%rah1
    real(rkx), pointer, contiguous, dimension(:,:) :: br => null( )        ! sfs%br
    real(rkx), pointer, contiguous, dimension(:,:,:) :: uxatm => null( )   ! atms%ubx3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vxatm => null( )   ! atms%vbx3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: udatm => null( )   ! atms%ubd3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vdatm => null( )   ! atms%vbd3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tatm => null( )    ! atms%tb3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: patm => null( )    ! atms%pb3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rhoatm => null( )  ! atms%rhob3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: patmf => null( )   ! atms%pf3d
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qxatm => null( ) ! atms%qx
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tkests => null( )  ! atms%tke
    real(rkx), pointer, contiguous, dimension(:,:,:) :: thatm => null( )   ! atms%th3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tvatm => null( )   ! atms%tv3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: za => null( )      ! atms%za
    real(rkx), pointer, contiguous, dimension(:,:,:) :: zq => null( )      ! atms%zq
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dzq => null( )     ! atms%dzq
    real(rkx), pointer, contiguous, dimension(:,:) :: rhox2d => null( )    ! atms%rhox2d
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chib => null( )  ! atms%chib3d
    real(rkx), pointer, contiguous, dimension(:,:,:) :: chifxuw => null( ) ! chifxuw
    real(rkx), pointer, contiguous, dimension(:,:,:) :: drydepv => null( ) ! drydepv
    real(rkx), pointer, contiguous, dimension(:,:,:) :: heatrt => null( )  ! radiation
    real(rkx), pointer, contiguous, dimension(:,:) :: uz0 => null( )   ! MYJ SF layer
    real(rkx), pointer, contiguous, dimension(:,:) :: vz0 => null( )   ! MYJ SF layer
    real(rkx), pointer, contiguous, dimension(:,:) :: thz0 => null( )  ! MYJ SF layer
    real(rkx), pointer, contiguous, dimension(:,:) :: qz0 => null( )   ! MYJ SF layer
  end type mod_2_pbl

  type pbl_2_mod
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tten => null( )       ! aten%t
    real(rkx), pointer, contiguous, dimension(:,:,:) :: uten => null( )       ! aten%u
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vten => null( )       ! aten%v
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qxten => null( )    ! aten%qx
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tketen => null( )     ! aten%tke
    real(rkx), pointer, contiguous, dimension(:,:,:) :: uxten => null( )      ! uxten%u
    real(rkx), pointer, contiguous, dimension(:,:,:) :: vxten => null( )      ! uxten%v
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chiten => null( )   ! chiten
    real(rkx), pointer, contiguous, dimension(:,:,:) :: remdrd => null( )     ! remdrd
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tkepbl => null( )     ! MYJ pbl
    real(rkx), pointer, contiguous, dimension(:,:) :: zpbl => null( )
    integer(ik4), pointer, contiguous, dimension(:,:) :: kpbl => null( )
  end type pbl_2_mod

  type mod_2_micro
    real(rkx), pointer, contiguous, dimension(:,:) :: xlat => null( )     ! mddom
    real(rkx), pointer, contiguous, dimension(:,:) :: psb => null( )      ! sfc
    real(rkx), pointer, contiguous, dimension(:,:) :: ps2 => null( )      ! from atms
    real(rkx), pointer, contiguous, dimension(:,:,:) :: phs => null( )    ! from atms
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pfs => null( )    ! from atms
    real(rkx), pointer, contiguous, dimension(:,:,:) :: delz => null( )   ! from atms
    real(rkx), pointer, contiguous, dimension(:,:,:) :: t => null( )      ! from atms
    real(rkx), pointer, contiguous, dimension(:,:,:) :: z => null( )      ! from atms
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rho => null( )    ! from atms
    real(rkx), pointer, contiguous, dimension(:,:,:) :: pverv => null( )  ! from atms
    real(rkx), pointer, contiguous, dimension(:,:,:) :: verv => null( )   ! from atms
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qxx => null( )  ! from atms
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qs => null( )     ! from atms
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rh => null( )     ! from atms
    real(rkx), pointer, contiguous, dimension(:,:,:) :: heatrt => null( ) ! radiation
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qdetr => null( )  ! conv. detr.
    real(rkx), pointer, contiguous, dimension(:,:,:) :: ccn => null( )    ! CCN from chem
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qvn => null( )    ! qxx(iqv)
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qcn => null( )    ! qxx(iqc)
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qin => null( )    ! qxx(iqi)
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qsn => null( )    ! qxx(iqs)
    real(rkx), pointer, contiguous, dimension(:,:,:) :: qrn => null( )    ! qxx(iqr)
    real(rkx), pointer, contiguous, dimension(:,:,:) :: cldf => null( )   ! cloud fraction
    integer(ik4), pointer, contiguous, dimension(:,:) :: ldmsk => null( ) ! mddom
    integer(ik4), pointer, contiguous, dimension(:,:) :: iveg => null( )  ! mddom
  end type mod_2_micro

  type micro_2_mod
    real(rkx), pointer, contiguous, dimension(:,:,:) :: fcc => null( )   ! Cloud cover
    real(rkx), pointer, contiguous, dimension(:,:) :: trrate => null( )  ! trrate
    real(rkx), pointer, contiguous, dimension(:,:) :: rainnc => null( )  ! sfc
    real(rkx), pointer, contiguous, dimension(:,:) :: lsmrnc => null( )  ! sfc
    real(rkx), pointer, contiguous, dimension(:,:) :: snownc => null( )  ! sfc
    real(rkx), pointer, contiguous, dimension(:,:) :: grplnc => null( )  ! sfc
    real(rkx), pointer, contiguous, dimension(:,:) :: hailnc => null( )  ! sfc
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rainls => null( ) ! Rain here
    real(rkx), pointer, contiguous, dimension(:,:,:) :: remrat => null( ) ! Rain here
    real(rkx), pointer, contiguous, dimension(:,:,:) :: rembc => null( )  ! Rain here
    real(rkx), pointer, contiguous, dimension(:,:,:) :: tten => null( )    ! tendency
    real(rkx), pointer, contiguous, dimension(:,:,:,:) :: qxten => null( ) ! tendency
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dia_qcr => null( ) ! diag for ccn
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dia_qcl => null( ) ! diag for ccn
    real(rkx), pointer, contiguous, dimension(:,:,:) :: dia_acr => null( ) ! diag for ccn
  end type micro_2_mod

  type nogtom_stats
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statssupw => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statssupc => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statserosw => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statserosc => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statsdetrw => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statsdetrc => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statsevapw => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statsevapc => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statscond1w => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statscond1c => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statsdepos => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statsmelt => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statsfrz => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statsrainev => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statssnowev => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statsautocvw => null( )
    real(rkx), public, pointer, contiguous, dimension(:,:,:) :: statsautocvc => null( )
  end type nogtom_stats

end module mod_regcm_types

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
