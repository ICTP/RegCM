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
    logical , pointer , dimension(:,:) :: gmask => null( )
    logical , pointer , dimension(:,:,:) :: sgmask => null( )
    logical , pointer , dimension(:,:) :: global_gmask => null( )
    logical , pointer , dimension(:,:,:) :: global_sgmask => null( )
    logical , pointer , dimension(:,:) :: global_out_sgmask => null( )
    integer(ik4) , pointer , dimension(:) :: linear_npoint_g => null( )
    integer(ik4) , pointer , dimension(:) :: linear_displ_g => null( )
    integer(ik4) , pointer , dimension(:) :: cartesian_npoint_g => null( )
    integer(ik4) , pointer , dimension(:) :: cartesian_displ_g => null( )
    integer(ik4) , pointer , dimension(:) :: linear_npoint_sg => null( )
    integer(ik4) , pointer , dimension(:) :: linear_displ_sg => null( )
    integer(ik4) , pointer , dimension(:) :: cartesian_npoint_sg => null( )
    integer(ik4) , pointer , dimension(:) :: cartesian_displ_sg => null( )
  end type masked_comm

  type model_area
    logical :: bandflag
    logical :: crmflag
    logical :: has_bdy
    logical :: has_bdyleft , has_bdyright , has_bdytop , has_bdybottom
    logical :: has_bdytopleft , has_bdytopright
    logical :: has_bdybottomleft , has_bdybottomright
    integer(ik4) , dimension(2) :: location
    integer(ik4) :: left , right , top , bottom
    integer(ik4) :: topleft , topright , bottomleft , bottomright
    integer(ik4) :: ibt1 , ibt2 , ibt3 , ibt4 , ibt6
    integer(ik4) :: ibb1 , ibb2 , ibb3 , ibb4 , ibb6
    integer(ik4) :: jbl1 , jbl2 , jbl3 , jbl4 , jbl6
    integer(ik4) :: jbr1 , jbr2 , jbr3 , jbr4 , jbr6
  end type model_area

  type domain
    real(rkx) , pointer , dimension(:,:) :: ht => null( )
    real(rkx) , pointer , dimension(:,:) :: lndcat => null( )
    real(rkx) , pointer , dimension(:,:) :: lndtex => null( )
    real(rkx) , pointer , dimension(:,:) :: mask => null( )
    real(rkx) , pointer , dimension(:,:) :: area => null( )
    real(rkx) , pointer , dimension(:,:) :: dlat => null( )
    real(rkx) , pointer , dimension(:,:) :: dlon => null( )
    real(rkx) , pointer , dimension(:,:) :: ulat => null( )
    real(rkx) , pointer , dimension(:,:) :: ulon => null( )
    real(rkx) , pointer , dimension(:,:) :: vlat => null( )
    real(rkx) , pointer , dimension(:,:) :: vlon => null( )
    real(rkx) , pointer , dimension(:,:) :: xlat => null( )
    real(rkx) , pointer , dimension(:,:) :: xlon => null( )
    real(rkx) , pointer , dimension(:,:) :: msfu => null( )
    real(rkx) , pointer , dimension(:,:) :: msfv => null( )
    real(rkx) , pointer , dimension(:,:) :: hx => null( )
    real(rkx) , pointer , dimension(:,:) :: hy => null( )
    real(rkx) , pointer , dimension(:,:) :: msfx => null( )
    real(rkx) , pointer , dimension(:,:) :: msfd => null( )
    real(rkx) , pointer , dimension(:,:) :: coriol => null( )
    real(rkx) , pointer , dimension(:,:) :: coriou => null( )
    real(rkx) , pointer , dimension(:,:) :: coriov => null( )
    real(rkx) , pointer , dimension(:,:) :: ef => null( )
    real(rkx) , pointer , dimension(:,:) :: ddx => null( )
    real(rkx) , pointer , dimension(:,:) :: ddy => null( )
    real(rkx) , pointer , dimension(:,:) :: ex => null( )
    real(rkx) , pointer , dimension(:,:) :: crx => null( )
    real(rkx) , pointer , dimension(:,:) :: cry => null( )
    real(rkx) , pointer , dimension(:,:) :: dmdy => null( )
    real(rkx) , pointer , dimension(:,:) :: dmdx => null( )
    real(rkx) , pointer , dimension(:,:) :: xmsf => null( )
    real(rkx) , pointer , dimension(:,:) :: dmsf => null( )
    real(rkx) , pointer , dimension(:,:) :: snowam => null( )
    real(rkx) , pointer , dimension(:,:) :: smoist => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rmoist => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rts => null( )
    real(rkx) , pointer , dimension(:,:) :: dhlake => null( )
    integer(ik4) , pointer , dimension(:,:) :: ldmsk => null( )
    integer(ik4) , pointer , dimension(:,:) :: iveg => null( )
    integer(ik4) , pointer , dimension(:,:) :: itex => null( )
  end type domain

  type domain_subgrid
    real(rkx) , pointer , dimension(:,:,:) :: ht => null( )
    real(rkx) , pointer , dimension(:,:,:) :: lndcat => null( )
    real(rkx) , pointer , dimension(:,:,:) :: lndtex => null( )
    real(rkx) , pointer , dimension(:,:,:) :: xlat => null( )
    real(rkx) , pointer , dimension(:,:,:) :: xlon => null( )
    real(rkx) , pointer , dimension(:,:,:) :: mask => null( )
    real(rkx) , pointer , dimension(:,:,:) :: area => null( )
    real(rkx) , pointer , dimension(:,:,:) :: dhlake => null( )
    integer(ik4) , pointer , dimension(:,:,:) :: ldmsk => null( )
    integer(ik4) , pointer , dimension(:,:,:) :: iveg => null( )
    integer(ik4) , pointer , dimension(:,:,:) :: itex => null( )
  end type domain_subgrid

  type mass_divergence
    real(rkx) , pointer , dimension(:,:,:) :: cr => null( ) ! cross points
  end type mass_divergence

  type atmosphere
    real(rkx) , pointer , dimension(:,:,:) :: u => null( )
    real(rkx) , pointer , dimension(:,:,:) :: v => null( )
    real(rkx) , pointer , dimension(:,:,:) :: ux => null( )
    real(rkx) , pointer , dimension(:,:,:) :: vx => null( )
    real(rkx) , pointer , dimension(:,:,:) :: w => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pai => null( )
    real(rkx) , pointer , dimension(:,:,:) :: t => null( )
    real(rkx) , pointer , dimension(:,:,:) :: p => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rho => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tvirt => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tetav => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tke => null( )
    real(rkx) , pointer , dimension(:,:,:) :: qs => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: qx => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: trac => null( )
    real(rkx) , pointer , dimension(:,:,:) :: zeta => null( )
    real(rkx) , pointer , dimension(:,:,:) :: zetaf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: dz => null( )
    real(rkx) , pointer , dimension(:,:,:) :: fmz => null( )
    real(rkx) , pointer , dimension(:,:,:) :: fmzf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tten => null( )
    real(rkx) , pointer , dimension(:,:,:) :: uten => null( )
    real(rkx) , pointer , dimension(:,:,:) :: vten => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tketen => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: qxten => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: chiten => null( )
  end type atmosphere

  type reference_atmosphere
    real(rkx) , pointer , dimension(:,:) :: ps => null( )
    real(rkx) , pointer , dimension(:,:) :: psdot => null( )
    real(rkx) , pointer , dimension(:,:,:) :: t => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pr => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rho => null( )
    real(rkx) , pointer , dimension(:,:,:) :: z => null( )
    real(rkx) , pointer , dimension(:,:,:) :: zd => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rhof => null( )
    real(rkx) , pointer , dimension(:,:,:) :: zf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: dzf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: dprddx => null( )
    real(rkx) , pointer , dimension(:,:,:) :: dprddy => null( )
  end type reference_atmosphere

  type atmstate_a
    real(rkx) , pointer , dimension(:,:,:) :: u => null( )
    real(rkx) , pointer , dimension(:,:,:) :: v => null( )
    real(rkx) , pointer , dimension(:,:,:) :: w => null( )
    real(rkx) , pointer , dimension(:,:,:) :: t => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tke => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pp => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pr => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rho => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: qx => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: chi => null( )
  end type atmstate_a

  type atmstate_b
    real(rkx) , pointer , dimension(:,:,:) :: u => null( )
    real(rkx) , pointer , dimension(:,:,:) :: v => null( )
    real(rkx) , pointer , dimension(:,:,:) :: w => null( )
    real(rkx) , pointer , dimension(:,:,:) :: t => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tke => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pp => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pr => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: qx => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: chi => null( )
  end type atmstate_b

  type atmstate_c
    real(rkx) , pointer , dimension(:,:,:) :: u => null( )
    real(rkx) , pointer , dimension(:,:,:) :: v => null( )
    real(rkx) , pointer , dimension(:,:,:) :: w => null( )
    real(rkx) , pointer , dimension(:,:,:) :: t => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tke => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pp => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: qx => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: chi => null( )
  end type atmstate_c

  type atmstate_tendency
    real(rkx) , pointer , dimension(:,:,:,:) :: u => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: v => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: w => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: t => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: tke => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: pp => null( )
    real(rkx) , pointer , dimension(:,:,:,:,:) :: qx => null( )
    real(rkx) , pointer , dimension(:,:,:,:,:) :: chi => null( )
  end type atmstate_tendency

  type crosswind_tendency
    real(rkx) , pointer , dimension(:,:,:) :: u => null( )
    real(rkx) , pointer , dimension(:,:,:) :: v => null( )
    real(rkx) , pointer , dimension(:,:,:) :: ud => null( )
    real(rkx) , pointer , dimension(:,:,:) :: vd => null( )
  end type crosswind_tendency

  type atmstate_decoupled
    real(rkx) , pointer , dimension(:,:,:) :: uc => null( )  ! pressure
    real(rkx) , pointer , dimension(:,:,:) :: vc => null( )  ! pressure
    real(rkx) , pointer , dimension(:,:,:) :: umc => null( ) ! mapfactor
    real(rkx) , pointer , dimension(:,:,:) :: vmc => null( ) ! mapfactor
    real(rkx) , pointer , dimension(:,:,:) :: ud => null( )  ! De-coupled
    real(rkx) , pointer , dimension(:,:,:) :: vd => null( )  ! De-coupled
    real(rkx) , pointer , dimension(:,:,:) :: umd => null( ) ! mapfactor
    real(rkx) , pointer , dimension(:,:,:) :: vmd => null( ) ! mapfactor
    real(rkx) , pointer , dimension(:,:,:) :: w => null( )
    real(rkx) , pointer , dimension(:,:,:) :: t => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tv => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: qx => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pp => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pr => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rho => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: chi => null( )
  end type atmstate_decoupled

  type tcm_state
    ! Vertical momentum diffusivity
    real(rkx) , pointer , dimension(:,:,:) :: kzm => null( )    ! (m^2/s)
    ! Vertical scalar diffusivity
    real(rkx) , pointer , dimension(:,:,:) :: kth => null( )    ! (m^2/s)
  end type tcm_state

  type tendiag
    real(rkx) , pointer , dimension(:,:,:) :: adh => null( )
    real(rkx) , pointer , dimension(:,:,:) :: adv => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tbl => null( )
    real(rkx) , pointer , dimension(:,:,:) :: dif => null( )
    real(rkx) , pointer , dimension(:,:,:) :: bdy => null( )
    real(rkx) , pointer , dimension(:,:,:) :: con => null( )
    real(rkx) , pointer , dimension(:,:,:) :: adi => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rad => null( )
    real(rkx) , pointer , dimension(:,:,:) :: lsc => null( )
  end type tendiag

  type qendiag
    real(rkx) , pointer , dimension(:,:,:) :: adh => null( )
    real(rkx) , pointer , dimension(:,:,:) :: adv => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tbl => null( )
    real(rkx) , pointer , dimension(:,:,:) :: dif => null( )
    real(rkx) , pointer , dimension(:,:,:) :: bdy => null( )
    real(rkx) , pointer , dimension(:,:,:) :: con => null( )
    real(rkx) , pointer , dimension(:,:,:) :: adi => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rad => null( )
    real(rkx) , pointer , dimension(:,:,:) :: lsc => null( )
    real(rkx) , pointer , dimension(:,:,:) :: qcl => null( )
    real(rkx) , pointer , dimension(:,:,:) :: qcr => null( )
    real(rkx) , pointer , dimension(:,:,:) :: acr => null( )
  end type qendiag

  type surfstate
    real(rkx) , pointer , dimension(:,:) :: psa => null( )
    real(rkx) , pointer , dimension(:,:) :: psb => null( )
    real(rkx) , pointer , dimension(:,:) :: psc => null( )
    real(rkx) , pointer , dimension(:,:) :: psdota => null( )
    real(rkx) , pointer , dimension(:,:) :: psdotb => null( )
    real(rkx) , pointer , dimension(:,:) :: tg => null( )
    real(rkx) , pointer , dimension(:,:) :: rainc => null( )
    real(rkx) , pointer , dimension(:,:) :: rainnc => null( )
    real(rkx) , pointer , dimension(:,:) :: snownc => null( )
    real(rkx) , pointer , dimension(:,:) :: grplnc => null( )
    real(rkx) , pointer , dimension(:,:) :: hailnc => null( )
    real(rkx) , pointer , dimension(:,:) :: hfx => null( )
    real(rkx) , pointer , dimension(:,:) :: qfx => null( )
    real(rkx) , pointer , dimension(:,:) :: tgbb => null( )
    real(rkx) , pointer , dimension(:,:) :: q2m => null( )
    real(rkx) , pointer , dimension(:,:) :: uvdrag => null( )
    real(rkx) , pointer , dimension(:,:) :: ustar => null( )
    real(rkx) , pointer , dimension(:,:) :: u10m => null( )
    real(rkx) , pointer , dimension(:,:) :: v10m => null( )
    real(rkx) , pointer , dimension(:,:) :: w10m => null( )
    real(rkx) , pointer , dimension(:,:) :: zo => null( )
    real(rkx) , pointer , dimension(:,:) :: ram1 => null( )
    real(rkx) , pointer , dimension(:,:) :: rah1 => null( )
    real(rkx) , pointer , dimension(:,:) :: br => null( )
    real(rkx) , pointer , dimension(:,:) :: uz0 => null( )    ! MYJ SF layer
    real(rkx) , pointer , dimension(:,:) :: vz0 => null( )    ! MYJ SF layer
    real(rkx) , pointer , dimension(:,:) :: thz0 => null( )   ! MYJ SF layer
    real(rkx) , pointer , dimension(:,:) :: qz0 => null( )    ! MYJ SF layer
    real(rkx) , pointer , dimension(:,:) :: dtrnof => null( ) ! rnoff
  end type surfstate

  type slice
    real(rkx) , pointer , dimension(:,:,:) :: th3d => null( )
    real(rkx) , pointer , dimension(:,:) :: th700 => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tb3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tv3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pb3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: pf3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rhob3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: ubx3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: vbx3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: wb3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: ppb3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: wpx3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: ubd3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: vbd3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rhb3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: qsb3d => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: qxb3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: zq => null( )
    real(rkx) , pointer , dimension(:,:,:) :: za => null( )
    real(rkx) , pointer , dimension(:,:,:) :: dzq => null( )
    real(rkx) , pointer , dimension(:,:) :: rhox2d => null( )
    real(rkx) , pointer , dimension(:,:) :: tp2d => null( )
    real(rkx) , pointer , dimension(:,:) :: ps2d => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: chib3d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tkepbl => null( )
  end type slice

  type v3dbound
    real(rkx) , pointer , dimension(:,:,:) :: b0 => null( )
    real(rkx) , pointer , dimension(:,:,:) :: b1 => null( )
    real(rkx) , pointer , dimension(:,:,:) :: bt => null( )
  end type v3dbound

  type v2dbound
    real(rkx) , pointer , dimension(:,:) :: b0 => null( )
    real(rkx) , pointer , dimension(:,:) :: b1 => null( )
    real(rkx) , pointer , dimension(:,:) :: bt => null( )
  end type v2dbound

  type nhboundhelp
    real(rkx) , pointer , dimension(:,:,:) :: tvirt => null( )
    real(rkx) , pointer , dimension(:,:) :: ps => null( )
  end type nhboundhelp

  type bound_area
    logical :: havebound
    logical , pointer , dimension(:,:) :: bsouth => null( )
    logical , pointer , dimension(:,:) :: bnorth => null( )
    logical , pointer , dimension(:,:) :: beast => null( )
    logical , pointer , dimension(:,:) :: bwest => null( )
    integer(ik4) :: ns , nn , ne , nw
    integer(ik4) :: nsp
    integer(ik4) , pointer , dimension(:,:) :: ibnd => null( )
  end type bound_area

  type exp_data
    real(rk8) , pointer , dimension(:,:) :: psfc => null( )
    real(rk8) , pointer , dimension(:,:) :: tsfc => null( )
    real(rk8) , pointer , dimension(:,:) :: qsfc => null( )
    real(rk8) , pointer , dimension(:,:) :: swrd => null( )
    real(rk8) , pointer , dimension(:,:) :: lwrd => null( )
    real(rk8) , pointer , dimension(:,:) :: dlwr => null( )
    real(rk8) , pointer , dimension(:,:) :: lhfx => null( )
    real(rk8) , pointer , dimension(:,:) :: shfx => null( )
    real(rk8) , pointer , dimension(:,:) :: prec => null( )
    real(rk8) , pointer , dimension(:,:) :: wndu => null( )
    real(rk8) , pointer , dimension(:,:) :: wndv => null( )
    real(rk8) , pointer , dimension(:,:) :: rnof => null( )
    real(rk8) , pointer , dimension(:,:) :: snof => null( )
    real(rk8) , pointer , dimension(:,:) :: taux => null( )
    real(rk8) , pointer , dimension(:,:) :: tauy => null( )
    real(rk8) , pointer , dimension(:,:) :: wspd => null( )
    real(rk8) , pointer , dimension(:,:) :: wdir => null( )
    real(rk8) , pointer , dimension(:,:) :: ustr => null( )
    real(rk8) , pointer , dimension(:,:) :: nflx => null( )
    real(rk8) , pointer , dimension(:,:) :: sflx => null( )
    real(rk8) , pointer , dimension(:,:) :: snow => null( )
    real(rk8) , pointer , dimension(:,:) :: dswr => null( )
    real(rk8) , pointer , dimension(:,:) :: rhoa => null( )
  end type exp_data

  type exp_data3d
    real(rk8) , pointer , dimension(:,:,:) :: t => null( )
    real(rk8) , pointer , dimension(:,:,:) :: q => null( )
    real(rk8) , pointer , dimension(:,:,:) :: u => null( )
    real(rk8) , pointer , dimension(:,:,:) :: v => null( )
    real(rk8) , pointer , dimension(:,:,:) :: w => null( )
    real(rk8) , pointer , dimension(:,:,:) :: cldfrc => null( )
    real(rk8) , pointer , dimension(:,:,:) :: cldlwc => null( )
  end type exp_data3d

  type imp_data
    real(rkx) , pointer , dimension(:,:) :: sst => null( )
    real(rkx) , pointer , dimension(:,:) :: sit => null( )
    real(rkx) , pointer , dimension(:,:) :: msk => null( )
    real(rkx) , pointer , dimension(:,:) :: zo => null( )
    real(rkx) , pointer , dimension(:,:) :: ustar => null( )
  end type imp_data

  type lm_state
    real(rkx) , pointer , dimension(:,:,:) :: gwet => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: sw => null( )
    real(rkx) , pointer , dimension(:,:,:) :: ssw => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rsw => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tsw => null( )
    real(rkx) , pointer , dimension(:,:,:) :: lncl => null( )
    real(rkx) , pointer , dimension(:,:,:) :: ldew => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tgbb => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tgrd => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tgbrd => null( )
    real(rkx) , pointer , dimension(:,:,:) :: taf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tlef => null( )
    real(rkx) , pointer , dimension(:,:,:) :: sfice => null( )
    real(rkx) , pointer , dimension(:,:,:) :: snag => null( )
    real(rkx) , pointer , dimension(:,:,:) :: sncv => null( )
    real(rkx) , pointer , dimension(:,:,:) :: scvk => null( )
    real(rkx) , pointer , dimension(:,:,:) :: emisv => null( )
    real(rkx) , pointer , dimension(:,:,:) :: sent => null( )
    real(rkx) , pointer , dimension(:,:,:) :: evpr => null( )
    real(rkx) , pointer , dimension(:,:,:) :: deltat => null( )
    real(rkx) , pointer , dimension(:,:,:) :: deltaq => null( )
    real(rkx) , pointer , dimension(:,:,:) :: drag => null( )
    real(rkx) , pointer , dimension(:,:,:) :: prcp => null( )
    real(rkx) , pointer , dimension(:,:,:) :: snwm => null( )
    real(rkx) , pointer , dimension(:,:,:) :: trnof => null( )
    real(rkx) , pointer , dimension(:,:,:) :: sigf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: sfcp => null( )
    real(rkx) , pointer , dimension(:,:,:) :: srnof => null( )
    real(rkx) , pointer , dimension(:,:,:) :: xlai => null( )
    real(rkx) , pointer , dimension(:,:,:) :: q2m => null( )
    real(rkx) , pointer , dimension(:,:,:) :: t2m => null( )
    real(rkx) , pointer , dimension(:,:,:) :: u10m => null( )
    real(rkx) , pointer , dimension(:,:,:) :: v10m => null( )
    real(rkx) , pointer , dimension(:,:,:) :: ram1 => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rah1 => null( )
    real(rkx) , pointer , dimension(:,:,:) :: br => null( )
    real(rkx) , pointer , dimension(:,:,:) :: taux => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tauy => null( )
    real(rkx) , pointer , dimension(:,:,:) :: swalb => null( )
    real(rkx) , pointer , dimension(:,:,:) :: lwalb => null( )
    real(rkx) , pointer , dimension(:,:,:) :: urlwf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: swdiralb => null( )
    real(rkx) , pointer , dimension(:,:,:) :: lwdiralb => null( )
    real(rkx) , pointer , dimension(:,:,:) :: swdifalb => null( )
    real(rkx) , pointer , dimension(:,:,:) :: lwdifalb => null( )
    real(rkx) , pointer , dimension(:,:,:) :: wt => null( )
    real(rkx) , pointer , dimension(:,:,:) :: eta => null( )
    real(rkx) , pointer , dimension(:,:,:) :: hi => null( )
    real(rkx) , pointer , dimension(:,:,:) :: hsnow => null( )
    real(rkx) , pointer , dimension(:,:,:) :: um10 => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: tlake => null( )
    logical , pointer , dimension(:,:,:) :: lakmsk => null( )
    real(rkx) , pointer , dimension(:,:,:) :: deltas => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tdeltas => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tskin => null( )
    real(rkx) , pointer , dimension(:,:,:) :: sst => null( )
    real(rkx) , pointer , dimension(:,:,:) :: zo => null( )
    real(rkx) , pointer , dimension(:,:,:) :: ustar => null( )
    real(rkx) , pointer , dimension(:,:,:) :: w10m => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rhoa => null( )
#ifdef CLM45
    real(rkx) , pointer , dimension(:,:,:) :: hfso => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: vocemiss => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: dustemiss => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: ddepv => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: sw_vol => null( )
    real(rkx) , pointer , dimension(:,:,:,:) :: tsoi => null( )
#endif
  end type lm_state

  type lm_exchange
    real(rkx) , pointer , dimension(:,:) :: ssw2da => null( )
    real(rkx) , pointer , dimension(:,:) :: sfracv2d => null( )
    real(rkx) , pointer , dimension(:,:) :: sfracb2d => null( )
    real(rkx) , pointer , dimension(:,:) :: sfracs2d => null( )
    real(rkx) , pointer , dimension(:,:) :: svegfrac2d => null( )
    real(rkx) , pointer , dimension(:,:) :: sxlai2d => null( )
    real(rkx) , pointer , dimension(:,:,:) :: sw_vol => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tsoi => null( )
    real(rkx) , pointer , dimension(:,:) :: xlat => null( )     ! mddom%xlat
    real(rkx) , pointer , dimension(:,:) :: xlon => null( )     ! mddom%xlon
    real(rkx) , pointer , dimension(:,:) :: lndcat => null( )   ! mddom%lndcat
    real(rkx) , pointer , dimension(:,:) :: ht => null( )       ! mddom%ht
    real(rkx) , pointer , dimension(:,:) :: snowam => null( )   ! mddom%snowam
    real(rkx) , pointer , dimension(:,:) :: smoist => null( )   ! mddom%smoist
    real(rkx) , pointer , dimension(:,:,:) :: rmoist => null( ) ! mddom%rmoist
    real(rkx) , pointer , dimension(:,:,:) :: rts => null( )    ! mddom%rts
    integer(ik4) , pointer , dimension(:,:) :: iveg => null( )  ! mddom%iveg
    integer(ik4) , pointer , dimension(:,:) :: itex => null( )  ! mddom%itex
    integer(ik4) , pointer , dimension(:,:) :: ldmsk => null( ) ! mddom%ldmsk
    real(rkx) , pointer , dimension(:,:,:) :: ht1 => null( )    ! mdsub%ht
    real(rkx) , pointer , dimension(:,:,:) :: lndcat1 => null( ) ! mdsub%lndcat
    real(rkx) , pointer , dimension(:,:,:) :: area1 => null( )   ! mdsub%area
    real(rkx) , pointer , dimension(:,:,:) :: xlat1 => null( )   ! mdsub%xlat
    real(rkx) , pointer , dimension(:,:,:) :: xlon1 => null( )   ! mdsub%xlon
    real(rkx) , pointer , dimension(:,:,:) :: dhlake1 => null( ) ! mdsub%dhlake
    integer(ik4) , pointer , dimension(:,:,:) :: ldmsk1 => null( ) ! mdsub%ldmsk
    integer(ik4) , pointer , dimension(:,:,:) :: iveg1 => null( )  ! mdsub%iveg
    integer(ik4) , pointer , dimension(:,:,:) :: itex1 => null( )  ! mdsub%itex
    integer(ik4) , pointer , dimension(:,:) :: icplmsk => null( )  ! cplmsk
    real(rkx) , pointer , dimension(:,:) :: patm => null( ) ! atms%pb3d(kz)
    real(rkx) , pointer , dimension(:,:) :: uatm => null( ) ! atms%ubx3d(kz)
    real(rkx) , pointer , dimension(:,:) :: vatm => null( ) ! atms%vbx3d(kz)
    real(rkx) , pointer , dimension(:,:) :: tatm => null( ) ! atms%tb3d(kz)
    real(rkx) , pointer , dimension(:,:) :: thatm => null( )! atms%th3d(kz)
    real(rkx) , pointer , dimension(:,:) :: qvatm => null( )! atms%qxb3d(kz,iqv)
    real(rkx) , pointer , dimension(:,:) :: hgt => null( )  ! za(kz)
    real(rkx) , pointer , dimension(:,:) :: hpbl => null( ) ! zpbl
    real(rkx) , pointer , dimension(:,:) :: hfx => null( )  ! sfs%hfx
    real(rkx) , pointer , dimension(:,:) :: qfx => null( )  ! sfs%qfx
    real(rkx) , pointer , dimension(:,:) :: totcf => null( ) ! totcf
    real(rkx) , pointer , dimension(:,:) :: sfps => null( )  ! sfs%psb
    real(rkx) , pointer , dimension(:,:) :: sfta => null( )  ! atms%ts2d
    real(rkx) , pointer , dimension(:,:) :: uvdrag => null( ) ! sfs%uvdrag
    real(rkx) , pointer , dimension(:,:) :: tg => null( )     ! sfs%tg
    real(rkx) , pointer , dimension(:,:) :: tgbb => null( )   ! sfs%tgbb
    real(rkx) , pointer , dimension(:,:) :: ram1 => null( )   ! sfs%ram1
    real(rkx) , pointer , dimension(:,:) :: rah1 => null( )   ! sfs%rah1
    real(rkx) , pointer , dimension(:,:) :: br => null( )     ! sfs%br
    real(rkx) , pointer , dimension(:,:) :: u10m => null( )   ! sfs%u10m
    real(rkx) , pointer , dimension(:,:) :: v10m => null( )   ! sfs%v10m
    real(rkx) , pointer , dimension(:,:) :: q2m => null( )    ! sfs%q2m
    real(rkx) , pointer , dimension(:,:) :: rhox => null( )   ! rhox2d
    real(rkx) , pointer , dimension(:,:) :: rswf => null( )   ! fsw
    real(rkx) , pointer , dimension(:,:) :: rlwf => null( )   ! flw
    real(rkx) , pointer , dimension(:,:) :: dwrlwf => null( ) ! flwd
    real(rkx) , pointer , dimension(:,:) :: zencos => null( )   ! coszrs
    real(rkx) , pointer , dimension(:,:) :: dtrnof => null( ) ! rnoff
    real(rkx) , pointer , dimension(:,:) :: ncprate => null( )  ! pptnc
    real(rkx) , pointer , dimension(:,:) :: cprate => null( )   ! cprate
    real(rkx) , pointer , dimension(:,:) :: snwrat => null( )   ! snwrat
    real(rkx) , pointer , dimension(:,:) :: csrate => null( )   ! csrate
    real(rkx) , pointer , dimension(:,:) :: grprat => null( )   ! grprat
    real(rkx) , pointer , dimension(:,:) :: hairat => null( )   ! hairat
    real(rkx) , pointer , dimension(:,:) :: vegswab => null( )  ! sabveg
    real(rkx) , pointer , dimension(:,:) :: lwalb => null( )    ! albvl
    real(rkx) , pointer , dimension(:,:) :: swalb => null( )    ! albvs
    real(rkx) , pointer , dimension(:,:) :: swdiralb => null( ) ! aldirs
    real(rkx) , pointer , dimension(:,:) :: swdifalb => null( ) ! aldifs
    real(rkx) , pointer , dimension(:,:) :: lwdiralb => null( ) ! aldirl
    real(rkx) , pointer , dimension(:,:) :: lwdifalb => null( ) ! aldifl
    real(rkx) , pointer , dimension(:,:) :: swdir => null( )    ! solvs
    real(rkx) , pointer , dimension(:,:) :: swdif => null( )    ! solvsd
    real(rkx) , pointer , dimension(:,:) :: lwdir => null( )    ! solvl
    real(rkx) , pointer , dimension(:,:) :: lwdif => null( )    ! solvld
    real(rkx) , pointer , dimension(:,:) :: solinc => null( )   ! sinc
    real(rkx) , pointer , dimension(:,:) :: solar => null( )    ! solis
    real(rkx) , pointer , dimension(:,:) :: emissivity => null( ) ! emiss
    real(rkx) , pointer , dimension(:,:) :: deltaq => null( )     ! sdelq
    real(rkx) , pointer , dimension(:,:) :: deltat => null( )     ! sdelt
    real(rkx) , pointer , dimension(:,:,:) :: drydepflx => null( ) ! drydepflx
    real(rkx) , pointer , dimension(:,:,:) :: wetdepflx => null( ) ! wetdepflx
    integer(ik4) , pointer , dimension(:) :: idust => null( ) ! dust indices
    real(rkx) , pointer , dimension(:,:) :: zo => null( )     ! zo
    real(rkx) , pointer , dimension(:,:) :: ustar => null( )  ! ustar
    real(rkx) , pointer , dimension(:,:) :: w10m => null( )   ! w10m
#ifdef CLM
    real(rkx) , pointer , dimension(:,:,:) :: dep_vels => null( )
#endif
  end type lm_exchange

  type mod_2_rad
    real(rkx) , pointer , dimension(:,:,:) :: cldfrc => null( )
    real(rkx) , pointer , dimension(:,:,:) :: cldlwc => null( )
    real(rkx) , pointer , dimension(:,:,:) :: tatms => null( )     ! atms%tb3d
    real(rkx) , pointer , dimension(:,:,:) :: rhatms => null( )    ! atms%rhb3d
    real(rkx) , pointer , dimension(:,:,:) :: rhoatms => null( )   ! atms%rhob3d
    real(rkx) , pointer , dimension(:,:,:) :: phatms => null( )    ! atms%pb3d
    real(rkx) , pointer , dimension(:,:,:) :: pfatms => null( )    ! atms%pf3d
    real(rkx) , pointer , dimension(:,:,:) :: za => null( )        ! atms%za
    real(rkx) , pointer , dimension(:,:,:) :: zq => null( )        ! atms%zq
    real(rkx) , pointer , dimension(:,:,:) :: deltaz => null( )    ! atms%dzq
    real(rkx) , pointer , dimension(:,:) :: psatms => null( )      ! atms%ps2d
    real(rkx) , pointer , dimension(:,:,:,:) :: qxatms => null( )  ! atms%qxb3d
    real(rkx) , pointer , dimension(:,:,:,:) :: chiatms => null( ) ! atms%chib3d
    real(rkx) , pointer , dimension(:,:) :: tg => null( )     ! sfs%tgbb
    real(rkx) , pointer , dimension(:,:) :: xlat => null( )   ! mddom%xlat
    real(rkx) , pointer , dimension(:,:) :: xlon => null( )   ! mddom%xlon
    real(rkx) , pointer , dimension(:,:) :: ht => null( )     ! mddom%ht
    real(rkx) , pointer , dimension(:,:) :: ptrop => null( )
    real(rkx) , pointer , dimension(:,:) :: coszrs => null( )
    real(rkx) , pointer , dimension(:,:) :: albvs => null( )
    real(rkx) , pointer , dimension(:,:) :: albvl => null( )
    real(rkx) , pointer , dimension(:,:) :: aldirs => null( )
    real(rkx) , pointer , dimension(:,:) :: aldifs => null( )
    real(rkx) , pointer , dimension(:,:) :: aldirl => null( )
    real(rkx) , pointer , dimension(:,:) :: aldifl => null( )
    real(rkx) , pointer , dimension(:,:) :: emiss => null( )
    real(rkx) , pointer , dimension(:,:) :: ps0 => null( )
    real(rkx) , pointer , dimension(:,:) :: bps0 => null( )
    real(rkx) , pointer , dimension(:,:,:) :: btv0 => null( )
    real(rkx) , pointer , dimension(:,:) :: bps1 => null( )
    real(rkx) , pointer , dimension(:,:,:) :: btv1 => null( )
    integer(ik4) , pointer , dimension(:,:) :: ldmsk => null( )
  end type mod_2_rad

  type rad_2_mod
    real(rkx) , pointer , dimension(:,:) :: solis => null( )
    real(rkx) , pointer , dimension(:,:) :: sinc => null( )
    real(rkx) , pointer , dimension(:,:) :: sabveg => null( )
    real(rkx) , pointer , dimension(:,:) :: sols => null( )
    real(rkx) , pointer , dimension(:,:) :: soll => null( )
    real(rkx) , pointer , dimension(:,:) :: solvs => null( )
    real(rkx) , pointer , dimension(:,:) :: solvsd => null( )
    real(rkx) , pointer , dimension(:,:) :: solvl => null( )
    real(rkx) , pointer , dimension(:,:) :: solvld => null( )
    real(rkx) , pointer , dimension(:,:) :: fsw => null( )
    real(rkx) , pointer , dimension(:,:) :: flw => null( )
    real(rkx) , pointer , dimension(:,:) :: flwd => null( )
    real(rkx) , pointer , dimension(:,:) :: totcf => null( )
    real(rkx) , pointer , dimension(:,:,:) :: heatrt => null( )
  end type rad_2_mod

  type mod_2_cum
    real(rkx) , pointer , dimension(:,:) :: ht => null( )        ! mddom%ht
    real(rkx) , pointer , dimension(:,:) :: psa => null( )       ! sfs%psa
    real(rkx) , pointer , dimension(:,:) :: psb => null( )       ! sfs%psb
    real(rkx) , pointer , dimension(:,:) :: psdotb => null( )    ! sfs%psdotb
    real(rkx) , pointer , dimension(:,:) :: psf => null( )       ! atms%ps2d
    real(rkx) , pointer , dimension(:,:,:) :: pas => null( )     ! atms%pb3d
    real(rkx) , pointer , dimension(:,:,:) :: pasf => null( )    ! atms%pf3d
    real(rkx) , pointer , dimension(:,:,:) :: zas => null( )     ! atms%za
    real(rkx) , pointer , dimension(:,:,:) :: tas => null( )     ! atms%tb3d
    real(rkx) , pointer , dimension(:,:,:) :: uas => null( )     ! atms%ubx3d
    real(rkx) , pointer , dimension(:,:,:) :: vas => null( )     ! atms%vbx3d
    real(rkx) , pointer , dimension(:,:,:) :: was => null( )     ! atms%wx3d
    real(rkx) , pointer , dimension(:,:,:) :: wpas => null( )    ! atms%wpx3d
    real(rkx) , pointer , dimension(:,:,:) :: qsas => null( )    ! atms%qsb3d
    real(rkx) , pointer , dimension(:,:,:) :: tkeas => null( )   ! atms%tke
    real(rkx) , pointer , dimension(:,:,:) :: rhoas => null( )   ! atms%rhob3d
    real(rkx) , pointer , dimension(:,:,:) :: zfs => null( )     ! atms%zq
    real(rkx) , pointer , dimension(:,:,:) :: dzq => null( )     ! atms%dzq
    real(rkx) , pointer , dimension(:,:,:) :: qq1 => null( )     ! atm1%q
    real(rkx) , pointer , dimension(:,:,:) :: qdot => null( )    ! qdot
    real(rkx) , pointer , dimension(:,:,:,:) :: qxas => null( )  ! atms%qxb3d
    real(rkx) , pointer , dimension(:,:,:,:) :: chias => null( ) ! atms%chib3d
    real(rkx) , pointer , dimension(:,:) :: qfx => null( )       ! sfs%qfx
    real(rkx) , pointer , dimension(:,:) :: hfx => null( )       ! sfs%hfx
    real(rkx) , pointer , dimension(:,:,:) :: tten => null( )     ! aten%t
    real(rkx) , pointer , dimension(:,:,:) :: uten => null( )     ! aten%u
    real(rkx) , pointer , dimension(:,:,:) :: vten => null( )     ! aten%v
    real(rkx) , pointer , dimension(:,:,:,:) :: qxten => null( )  ! aten%qx
    real(rkx) , pointer , dimension(:,:,:,:) :: dynqx => null( )  ! aten%qx
    real(rkx) , pointer , dimension(:,:,:,:) :: chiten => null( ) ! aten%chi
    real(rkx) , pointer , dimension(:,:,:) :: heatrt => null( )   ! radiation
    real(rkx) , pointer , dimension(:,:,:) :: ccn => null( )      ! ccn
    integer(ik4) , pointer , dimension(:,:) :: ktrop => null( )
    integer(ik4) , pointer , dimension(:,:) :: ldmsk => null( )   ! mddom
  end type mod_2_cum

  type cum_2_mod
    real(rkx) , pointer , dimension(:,:,:) :: tten => null( )     ! aten%t
    real(rkx) , pointer , dimension(:,:,:) :: uten => null( )     ! aten%u
    real(rkx) , pointer , dimension(:,:,:) :: vten => null( )     ! aten%v
    real(rkx) , pointer , dimension(:,:,:,:) :: qxten => null( )  ! aten%qx
    real(rkx) , pointer , dimension(:,:,:,:) :: chiten => null( ) ! aten%chi
    real(rkx) , pointer , dimension(:,:) :: rainc => null( )
    real(rkx) , pointer , dimension(:,:) :: pcratec => null( )
    real(rkx) , pointer , dimension(:,:) :: sratec => null( )
    real(rkx) , pointer , dimension(:,:) :: trrate => null( )     ! trrate
    real(rkx) , pointer , dimension(:,:,:) :: convpr => null( )
    real(rkx) , pointer , dimension(:,:,:) :: cldfrc => null( )
    real(rkx) , pointer , dimension(:,:,:) :: cldlwc => null( )
    real(rkx) , pointer , dimension(:,:,:) :: q_detr => null( )
    real(rkx) , pointer , dimension(:,:,:) :: rain_cc => null( )
    integer(ik4) , pointer , dimension(:,:) :: kcumtop => null( )
    integer(ik4) , pointer , dimension(:,:) :: kcumbot => null( )
  end type cum_2_mod

  type mod_2_pbl
    real(rkx) , pointer , dimension(:,:) :: coriol => null( )    ! mddom%coriol
    integer(ik4) , pointer , dimension(:,:) :: ldmsk => null( )  ! mddom%ldmsk
    real(rkx) , pointer , dimension(:,:) :: ht => null( )        ! mddom%ht
    real(rkx) , pointer , dimension(:,:) :: psb => null( )       ! sfs%psb
    real(rkx) , pointer , dimension(:,:) :: psdotb => null( )    ! sfs%psdotb
    real(rkx) , pointer , dimension(:,:) :: tg => null( )        ! sfs%tgbb
    real(rkx) , pointer , dimension(:,:) :: q2m => null( )       ! sfs%q2m
    real(rkx) , pointer , dimension(:,:) :: u10m => null( )      ! sfs%u10m
    real(rkx) , pointer , dimension(:,:) :: v10m => null( )      ! sfs%v10m
    real(rkx) , pointer , dimension(:,:) :: qfx => null( )       ! sfs%qfx
    real(rkx) , pointer , dimension(:,:) :: hfx => null( )       ! sfs%hfx
    real(rkx) , pointer , dimension(:,:) :: uvdrag => null( )    ! sfs%uvdrag
    real(rkx) , pointer , dimension(:,:) :: zo => null( )        ! sfs%zo
    real(rkx) , pointer , dimension(:,:) :: ustar => null( )     ! sfs%ustar
    real(rkx) , pointer , dimension(:,:) :: ram1 => null( )      ! sfs%ram1
    real(rkx) , pointer , dimension(:,:) :: rah1 => null( )      ! sfs%rah1
    real(rkx) , pointer , dimension(:,:) :: br => null( )        ! sfs%br
    real(rkx) , pointer , dimension(:,:,:) :: uxatm => null( )   ! atms%ubx3d
    real(rkx) , pointer , dimension(:,:,:) :: vxatm => null( )   ! atms%vbx3d
    real(rkx) , pointer , dimension(:,:,:) :: udatm => null( )   ! atms%ubd3d
    real(rkx) , pointer , dimension(:,:,:) :: vdatm => null( )   ! atms%vbd3d
    real(rkx) , pointer , dimension(:,:,:) :: tatm => null( )    ! atms%tb3d
    real(rkx) , pointer , dimension(:,:,:) :: patm => null( )    ! atms%pb3d
    real(rkx) , pointer , dimension(:,:,:) :: rhoatm => null( )  ! atms%rhob3d
    real(rkx) , pointer , dimension(:,:,:) :: patmf => null( )   ! atms%pf3d
    real(rkx) , pointer , dimension(:,:,:,:) :: qxatm => null( ) ! atms%qx
    real(rkx) , pointer , dimension(:,:,:) :: tkests => null( )  ! atms%tke
    real(rkx) , pointer , dimension(:,:,:) :: thatm => null( )   ! atms%th3d
    real(rkx) , pointer , dimension(:,:,:) :: tvatm => null( )   ! atms%tv3d
    real(rkx) , pointer , dimension(:,:,:) :: za => null( )      ! atms%za
    real(rkx) , pointer , dimension(:,:,:) :: zq => null( )      ! atms%zq
    real(rkx) , pointer , dimension(:,:,:) :: dzq => null( )     ! atms%dzq
    real(rkx) , pointer , dimension(:,:) :: rhox2d => null( )    ! atms%rhox2d
    real(rkx) , pointer , dimension(:,:,:,:) :: chib => null( )  ! atms%chib3d
    real(rkx) , pointer , dimension(:,:,:) :: chifxuw => null( ) ! chifxuw
    real(rkx) , pointer , dimension(:,:,:) :: drydepv => null( ) ! drydepv
    real(rkx) , pointer , dimension(:,:,:) :: heatrt => null( )  ! radiation
    real(rkx) , pointer , dimension(:,:) :: uz0 => null( )   ! MYJ SF layer
    real(rkx) , pointer , dimension(:,:) :: vz0 => null( )   ! MYJ SF layer
    real(rkx) , pointer , dimension(:,:) :: thz0 => null( )  ! MYJ SF layer
    real(rkx) , pointer , dimension(:,:) :: qz0 => null( )   ! MYJ SF layer
  end type mod_2_pbl

  type pbl_2_mod
    real(rkx) , pointer , dimension(:,:,:) :: tten => null( )       ! aten%t
    real(rkx) , pointer , dimension(:,:,:) :: uten => null( )       ! aten%u
    real(rkx) , pointer , dimension(:,:,:) :: vten => null( )       ! aten%v
    real(rkx) , pointer , dimension(:,:,:,:) :: qxten => null( )    ! aten%qx
    real(rkx) , pointer , dimension(:,:,:) :: tketen => null( )     ! aten%tke
    real(rkx) , pointer , dimension(:,:,:) :: uxten => null( )      ! uxten%u
    real(rkx) , pointer , dimension(:,:,:) :: vxten => null( )      ! uxten%v
    real(rkx) , pointer , dimension(:,:,:,:) :: chiten => null( )   ! chiten
    real(rkx) , pointer , dimension(:,:,:) :: remdrd => null( )     ! remdrd
    real(rkx) , pointer , dimension(:,:,:) :: tkepbl => null( )     ! MYJ pbl
    real(rkx) , pointer , dimension(:,:) :: zpbl => null( )
    integer(ik4) , pointer , dimension(:,:) :: kpbl => null( )
  end type pbl_2_mod

  type mod_2_micro
    real(rkx) , pointer , dimension(:,:) :: xlat => null( )     ! mddom
    real(rkx) , pointer , dimension(:,:) :: psb => null( )      ! sfc
    real(rkx) , pointer , dimension(:,:) :: ps2 => null( )      ! from atms
    real(rkx) , pointer , dimension(:,:,:) :: phs => null( )    ! from atms
    real(rkx) , pointer , dimension(:,:,:) :: pfs => null( )    ! from atms
    real(rkx) , pointer , dimension(:,:,:) :: delz => null( )   ! from atms
    real(rkx) , pointer , dimension(:,:,:) :: t => null( )      ! from atms
    real(rkx) , pointer , dimension(:,:,:) :: z => null( )      ! from atms
    real(rkx) , pointer , dimension(:,:,:) :: rho => null( )    ! from atms
    real(rkx) , pointer , dimension(:,:,:) :: pverv => null( )  ! from atms
    real(rkx) , pointer , dimension(:,:,:) :: verv => null( )   ! from atms
    real(rkx) , pointer , dimension(:,:,:,:) :: qxx => null( )  ! from atms
    real(rkx) , pointer , dimension(:,:,:) :: qs => null( )     ! from atms
    real(rkx) , pointer , dimension(:,:,:) :: rh => null( )     ! from atms
    real(rkx) , pointer , dimension(:,:,:) :: heatrt => null( ) ! radiation
    real(rkx) , pointer , dimension(:,:,:) :: qdetr => null( )  ! conv. detr.
    real(rkx) , pointer , dimension(:,:,:) :: ccn => null( )    ! CCN from chem
    real(rkx) , pointer , dimension(:,:,:) :: qvn => null( )    ! qxx(iqv)
    real(rkx) , pointer , dimension(:,:,:) :: qcn => null( )    ! qxx(iqc)
    real(rkx) , pointer , dimension(:,:,:) :: qin => null( )    ! qxx(iqi)
    real(rkx) , pointer , dimension(:,:,:) :: qsn => null( )    ! qxx(iqs)
    real(rkx) , pointer , dimension(:,:,:) :: qrn => null( )    ! qxx(iqr)
    real(rkx) , pointer , dimension(:,:,:) :: cldf => null( )   ! cloud fraction
    integer(ik4) , pointer , dimension(:,:) :: ldmsk => null( ) ! mddom
    integer(ik4) , pointer , dimension(:,:) :: iveg => null( )  ! mddom
  end type mod_2_micro

  type micro_2_mod
    real(rkx) , pointer , dimension(:,:,:) :: fcc => null( )   ! Cloud cover
    real(rkx) , pointer , dimension(:,:) :: trrate => null( )  ! trrate
    real(rkx) , pointer , dimension(:,:) :: rainnc => null( )  ! sfc
    real(rkx) , pointer , dimension(:,:) :: lsmrnc => null( )  ! sfc
    real(rkx) , pointer , dimension(:,:) :: snownc => null( )  ! sfc
    real(rkx) , pointer , dimension(:,:) :: grplnc => null( )  ! sfc
    real(rkx) , pointer , dimension(:,:) :: hailnc => null( )  ! sfc
    real(rkx) , pointer , dimension(:,:,:) :: rainls => null( ) ! Rain here
    real(rkx) , pointer , dimension(:,:,:) :: remrat => null( ) ! Rain here
    real(rkx) , pointer , dimension(:,:,:) :: rembc => null( )  ! Rain here
    real(rkx) , pointer , dimension(:,:,:) :: tten => null( )    ! tendency
    real(rkx) , pointer , dimension(:,:,:,:) :: qxten => null( ) ! tendency
    real(rkx) , pointer , dimension(:,:,:) :: dia_qcr => null( ) ! diag for ccn
    real(rkx) , pointer , dimension(:,:,:) :: dia_qcl => null( ) ! diag for ccn
    real(rkx) , pointer , dimension(:,:,:) :: dia_acr => null( ) ! diag for ccn
  end type micro_2_mod

  type nogtom_stats
    real(rkx) , public  , pointer, dimension(:,:,:) :: statssupw => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statssupc => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statserosw => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statserosc => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statsdetrw => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statsdetrc => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statsevapw => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statsevapc => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statscond1w => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statscond1c => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statsdepos => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statsmelt => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statsfrz => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statsrainev => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statssnowev => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statsautocvw => null( )
    real(rkx) , public  , pointer, dimension(:,:,:) :: statsautocvc => null( )
  end type nogtom_stats

end module mod_regcm_types

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
