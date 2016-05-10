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
!
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
    logical , pointer , dimension(:,:) :: gmask
    logical , pointer , dimension(:,:,:) :: sgmask
    logical , pointer , dimension(:,:) :: global_gmask
    logical , pointer , dimension(:,:,:) :: global_sgmask
    logical , pointer , dimension(:,:) :: global_out_sgmask
    integer(ik4) , pointer , dimension(:) :: linear_npoint_g
    integer(ik4) , pointer , dimension(:) :: linear_displ_g
    integer(ik4) , pointer , dimension(:) :: cartesian_npoint_g
    integer(ik4) , pointer , dimension(:) :: cartesian_displ_g
    integer(ik4) , pointer , dimension(:) :: linear_npoint_sg
    integer(ik4) , pointer , dimension(:) :: linear_displ_sg
    integer(ik4) , pointer , dimension(:) :: cartesian_npoint_sg
    integer(ik4) , pointer , dimension(:) :: cartesian_displ_sg
  end type masked_comm

  type model_area
    logical :: bandflag
    logical :: has_bdy
    logical :: has_bdyleft , has_bdyright , has_bdytop , has_bdybottom
    logical :: has_bdytopleft , has_bdytopright
    logical :: has_bdybottomleft , has_bdybottomright
    integer(ik4) , dimension(2) :: location
    integer(ik4) :: left , right , top , bottom
    integer(ik4) :: topleft , topright , bottomleft , bottomright
    integer(ik4) :: ibt1 , ibt2 , ibt4 , ibt6 , ibb1 , ibb2 , ibb4 , ibb6
    integer(ik4) :: jbl1 , jbl2 , jbl4 , jbl6 , jbr1 , jbr2 , jbr4 , jbr6
  end type model_area

  type domain
    real(rkx) , pointer , dimension(:,:) :: ht
    real(rkx) , pointer , dimension(:,:) :: lndcat
    real(rkx) , pointer , dimension(:,:) :: xlat
    real(rkx) , pointer , dimension(:,:) :: xlon
    real(rkx) , pointer , dimension(:,:) :: mask
    real(rkx) , pointer , dimension(:,:) :: dlat
    real(rkx) , pointer , dimension(:,:) :: dlon
    real(rkx) , pointer , dimension(:,:) :: msfx
    real(rkx) , pointer , dimension(:,:) :: msfd
    real(rkx) , pointer , dimension(:,:) :: coriol
    real(rkx) , pointer , dimension(:,:) :: ef
    real(rkx) , pointer , dimension(:,:) :: ddx
    real(rkx) , pointer , dimension(:,:) :: ddy
    real(rkx) , pointer , dimension(:,:) :: ex
    real(rkx) , pointer , dimension(:,:) :: crx
    real(rkx) , pointer , dimension(:,:) :: cry
    real(rkx) , pointer , dimension(:,:) :: dmdy
    real(rkx) , pointer , dimension(:,:) :: dmdx
    real(rkx) , pointer , dimension(:,:) :: xmsf
    real(rkx) , pointer , dimension(:,:) :: dmsf
    real(rkx) , pointer , dimension(:,:) :: snowam
    real(rkx) , pointer , dimension(:,:) :: smoist
    real(rkx) , pointer , dimension(:,:,:) :: rmoist
    real(rkx) , pointer , dimension(:,:) :: dhlake
    integer(ik4) , pointer , dimension(:,:) :: ldmsk
    integer(ik4) , pointer , dimension(:,:) :: iveg
  end type domain

  type domain_subgrid
    real(rkx) , pointer , dimension(:,:,:) :: ht
    real(rkx) , pointer , dimension(:,:,:) :: lndcat
    real(rkx) , pointer , dimension(:,:,:) :: xlat
    real(rkx) , pointer , dimension(:,:,:) :: xlon
    real(rkx) , pointer , dimension(:,:,:) :: mask
    real(rkx) , pointer , dimension(:,:,:) :: dhlake
    integer(ik4) , pointer , dimension(:,:,:) :: ldmsk
    integer(ik4) , pointer , dimension(:,:,:) :: iveg
  end type domain_subgrid

  type mass_divergence
    real(rkx) , pointer , dimension(:,:,:) :: cr ! cross points
  end type mass_divergence

  type reference_atmosphere
    real(rkx) , pointer , dimension(:,:) :: ps
    real(rkx) , pointer , dimension(:,:) :: psdot
    real(rkx) , pointer , dimension(:,:,:) :: t
    real(rkx) , pointer , dimension(:,:,:) :: pr
    real(rkx) , pointer , dimension(:,:,:) :: pf
    real(rkx) , pointer , dimension(:,:,:) :: rho
  end type reference_atmosphere

  type atmstate_a
    real(rkx) , pointer , dimension(:,:,:) :: u
    real(rkx) , pointer , dimension(:,:,:) :: v
    real(rkx) , pointer , dimension(:,:,:) :: w
    real(rkx) , pointer , dimension(:,:,:) :: t
    real(rkx) , pointer , dimension(:,:,:,:) :: qx
    real(rkx) , pointer , dimension(:,:,:) :: tke
    real(rkx) , pointer , dimension(:,:,:) :: pp
    real(rkx) , pointer , dimension(:,:,:) :: pr
    real(rkx) , pointer , dimension(:,:,:) :: rho
  end type atmstate_a

  type atmstate_b
    real(rkx) , pointer , dimension(:,:,:) :: u
    real(rkx) , pointer , dimension(:,:,:) :: v
    real(rkx) , pointer , dimension(:,:,:) :: w
    real(rkx) , pointer , dimension(:,:,:) :: t
    real(rkx) , pointer , dimension(:,:,:,:) :: qx
    real(rkx) , pointer , dimension(:,:,:) :: tke
    real(rkx) , pointer , dimension(:,:,:) :: pp
    real(rkx) , pointer , dimension(:,:,:) :: pr
  end type atmstate_b

  type atmstate_c
    real(rkx) , pointer , dimension(:,:,:) :: u
    real(rkx) , pointer , dimension(:,:,:) :: v
    real(rkx) , pointer , dimension(:,:,:) :: w
    real(rkx) , pointer , dimension(:,:,:) :: t
    real(rkx) , pointer , dimension(:,:,:,:) :: qx
    real(rkx) , pointer , dimension(:,:,:) :: tke
    real(rkx) , pointer , dimension(:,:,:) :: pp
    real(rkx) , pointer , dimension(:,:,:) :: rho
  end type atmstate_c

  type atmstate_tendency
    real(rkx) , pointer , dimension(:,:,:) :: u
    real(rkx) , pointer , dimension(:,:,:) :: v
    real(rkx) , pointer , dimension(:,:,:) :: w
    real(rkx) , pointer , dimension(:,:,:) :: t
    real(rkx) , pointer , dimension(:,:,:,:) :: qx
    real(rkx) , pointer , dimension(:,:,:) :: tke
    real(rkx) , pointer , dimension(:,:,:) :: pp
  end type atmstate_tendency

  type atmstate_decoupled
    real(rkx) , pointer , dimension(:,:,:) :: uc  ! Coupled with pressure
    real(rkx) , pointer , dimension(:,:,:) :: vc  ! Coupled with pressure
    real(rkx) , pointer , dimension(:,:,:) :: umc ! Coupled * mapfactor
    real(rkx) , pointer , dimension(:,:,:) :: vmc ! Coupled * mapfactor
    real(rkx) , pointer , dimension(:,:,:) :: ud  ! De-coupled
    real(rkx) , pointer , dimension(:,:,:) :: vd  ! De-coupled
    real(rkx) , pointer , dimension(:,:,:) :: umd ! De-coupled * mapfactor
    real(rkx) , pointer , dimension(:,:,:) :: vmd ! De-coupled * mapfactor
    real(rkx) , pointer , dimension(:,:,:) :: w
    real(rkx) , pointer , dimension(:,:,:) :: t
    real(rkx) , pointer , dimension(:,:,:) :: tv
    real(rkx) , pointer , dimension(:,:,:,:) :: qx
    real(rkx) , pointer , dimension(:,:,:) :: pp
    real(rkx) , pointer , dimension(:,:,:) :: pr
    real(rkx) , pointer , dimension(:,:,:) :: rho
  end type atmstate_decoupled

  type tcm_state
    ! TKE*ps
    real(rkx) , pointer , dimension(:,:,:) :: tkeps  ! (m^2/s^2 * cb)
    ! Coupled TKE Advective Tendency
    real(rkx) , pointer , dimension(:,:,:) :: advtke ! (m^2/s^3 * cb)
    ! Vertical momentum diffusivity
    real(rkx) , pointer , dimension(:,:,:) :: kzm    ! (m^2/s)
    ! Vertical scalar diffusivity
    real(rkx) , pointer , dimension(:,:,:) :: kth    ! (m^2/s)
    ! Boundary layer height (m)
    real(rkx) , pointer , dimension(:,:) :: zpbl     ! (m)
    ! Surface layer TKE
    real(rkx) , pointer , dimension(:,:) :: srftke   ! (m^2/s^2)
  end type tcm_state

  type tendiag
    real(rkx) , pointer , dimension(:,:,:) :: adh
    real(rkx) , pointer , dimension(:,:,:) :: adv
    real(rkx) , pointer , dimension(:,:,:) :: tbl
    real(rkx) , pointer , dimension(:,:,:) :: dif
    real(rkx) , pointer , dimension(:,:,:) :: bdy
    real(rkx) , pointer , dimension(:,:,:) :: con
    real(rkx) , pointer , dimension(:,:,:) :: adi
    real(rkx) , pointer , dimension(:,:,:) :: rad
    real(rkx) , pointer , dimension(:,:,:) :: lsc
  end type tendiag

  type qendiag
    real(rkx) , pointer , dimension(:,:,:) :: adh
    real(rkx) , pointer , dimension(:,:,:) :: adv
    real(rkx) , pointer , dimension(:,:,:) :: tbl
    real(rkx) , pointer , dimension(:,:,:) :: dif
    real(rkx) , pointer , dimension(:,:,:) :: bdy
    real(rkx) , pointer , dimension(:,:,:) :: con
    real(rkx) , pointer , dimension(:,:,:) :: adi
    real(rkx) , pointer , dimension(:,:,:) :: rad
    real(rkx) , pointer , dimension(:,:,:) :: lsc
    real(rkx) , pointer , dimension(:,:,:) :: qcl
    real(rkx) , pointer , dimension(:,:,:) :: qcr
    real(rkx) , pointer , dimension(:,:,:) :: acr
  end type qendiag

  type surfstate
    real(rkx) , pointer , dimension(:,:) :: psa
    real(rkx) , pointer , dimension(:,:) :: psb
    real(rkx) , pointer , dimension(:,:) :: psc
    real(rkx) , pointer , dimension(:,:) :: psdota
    real(rkx) , pointer , dimension(:,:) :: psdotb
    real(rkx) , pointer , dimension(:,:) :: tga
    real(rkx) , pointer , dimension(:,:) :: tgb
    real(rkx) , pointer , dimension(:,:) :: rainc
    real(rkx) , pointer , dimension(:,:) :: rainnc
    real(rkx) , pointer , dimension(:,:) :: snownc
    real(rkx) , pointer , dimension(:,:) :: hfx
    real(rkx) , pointer , dimension(:,:) :: qfx
    real(rkx) , pointer , dimension(:,:) :: tgbb
    real(rkx) , pointer , dimension(:,:) :: uvdrag
    real(rkx) , pointer , dimension(:,:) :: ustar
    real(rkx) , pointer , dimension(:,:) :: zo
    real(rkx) , pointer , dimension(:,:) :: rhoa
  end type surfstate

  type slice
    real(rkx) , pointer , dimension(:,:,:) :: th3d
    real(rkx) , pointer , dimension(:,:) :: th700
    real(rkx) , pointer , dimension(:,:,:) :: tb3d
    real(rkx) , pointer , dimension(:,:,:) :: tp3d
    real(rkx) , pointer , dimension(:,:,:) :: pb3d
    real(rkx) , pointer , dimension(:,:,:) :: pf3d
    real(rkx) , pointer , dimension(:,:,:) :: rhob3d
    real(rkx) , pointer , dimension(:,:,:) :: ubx3d
    real(rkx) , pointer , dimension(:,:,:) :: vbx3d
    real(rkx) , pointer , dimension(:,:,:) :: wb3d
    real(rkx) , pointer , dimension(:,:,:) :: ppb3d
    real(rkx) , pointer , dimension(:,:,:) :: wpx3d
    real(rkx) , pointer , dimension(:,:,:) :: ubd3d
    real(rkx) , pointer , dimension(:,:,:) :: vbd3d
    real(rkx) , pointer , dimension(:,:,:) :: rhb3d
    real(rkx) , pointer , dimension(:,:,:) :: qsb3d
    real(rkx) , pointer , dimension(:,:,:,:) :: qxb3d
    real(rkx) , pointer , dimension(:,:,:) :: zq
    real(rkx) , pointer , dimension(:,:,:) :: za
    real(rkx) , pointer , dimension(:,:,:) :: dzq
    real(rkx) , pointer , dimension(:,:) :: rhox2d
    real(rkx) , pointer , dimension(:,:) :: ps2d
    real(rkx) , pointer , dimension(:,:,:,:) :: chib3d
  end type slice

  type diffx
    real(rkx) , pointer , dimension(:,:,:) :: t
    real(rkx) , pointer , dimension(:,:,:) :: u
    real(rkx) , pointer , dimension(:,:,:) :: v
    real(rkx) , pointer , dimension(:,:,:) :: w
    real(rkx) , pointer , dimension(:,:,:) :: pp
    real(rkx) , pointer , dimension(:,:,:,:) :: qx
  end type diffx

  type v3dbound
    real(rkx) , pointer , dimension(:,:,:) :: b0
    real(rkx) , pointer , dimension(:,:,:) :: b1
    real(rkx) , pointer , dimension(:,:,:) :: bt
  end type v3dbound

  type v2dbound
    real(rkx) , pointer , dimension(:,:) :: b0
    real(rkx) , pointer , dimension(:,:) :: b1
    real(rkx) , pointer , dimension(:,:) :: bt
  end type v2dbound

  type bound_area
    logical :: dotflag
    logical :: havebound
    logical , pointer , dimension(:,:) :: bsouth
    logical , pointer , dimension(:,:) :: bnorth
    logical , pointer , dimension(:,:) :: beast
    logical , pointer , dimension(:,:) :: bwest
    integer(ik4) :: ns , nn , ne , nw
    integer(ik4) :: nsp
    integer(ik4) , pointer , dimension(:,:) :: ibnd
  end type bound_area

  type exp_data
    real(rkx) , pointer , dimension(:,:) :: psfc
    real(rkx) , pointer , dimension(:,:) :: tsfc
    real(rkx) , pointer , dimension(:,:) :: qsfc
    real(rkx) , pointer , dimension(:,:) :: swrd
    real(rkx) , pointer , dimension(:,:) :: lwrd
    real(rkx) , pointer , dimension(:,:) :: dlwr
    real(rkx) , pointer , dimension(:,:) :: lhfx
    real(rkx) , pointer , dimension(:,:) :: shfx
    real(rkx) , pointer , dimension(:,:) :: prec
    real(rkx) , pointer , dimension(:,:) :: wndu
    real(rkx) , pointer , dimension(:,:) :: wndv
    real(rkx) , pointer , dimension(:,:) :: rnof
    real(rkx) , pointer , dimension(:,:) :: snof
    real(rkx) , pointer , dimension(:,:) :: taux
    real(rkx) , pointer , dimension(:,:) :: tauy
    real(rkx) , pointer , dimension(:,:) :: wspd
    real(rkx) , pointer , dimension(:,:) :: wdir
    real(rkx) , pointer , dimension(:,:) :: ustr
    real(rkx) , pointer , dimension(:,:) :: nflx
    real(rkx) , pointer , dimension(:,:) :: sflx
    real(rkx) , pointer , dimension(:,:) :: snow
    real(rkx) , pointer , dimension(:,:) :: dswr
    real(rkx) , pointer , dimension(:,:) :: rhoa
  end type exp_data

  type imp_data
    real(rkx) , pointer , dimension(:,:) :: sst
    real(rkx) , pointer , dimension(:,:) :: sit
    real(rkx) , pointer , dimension(:,:) :: msk
    real(rkx) , pointer , dimension(:,:) :: zo
    real(rkx) , pointer , dimension(:,:) :: ustar
  end type imp_data

  type lm_state
    real(rkx) , pointer , dimension(:,:,:) :: gwet
    real(rkx) , pointer , dimension(:,:,:,:) :: sw
    real(rkx) , pointer , dimension(:,:,:) :: ssw
    real(rkx) , pointer , dimension(:,:,:) :: rsw
    real(rkx) , pointer , dimension(:,:,:) :: tsw
    real(rkx) , pointer , dimension(:,:,:) :: ldew
    real(rkx) , pointer , dimension(:,:,:) :: lncl
    real(rkx) , pointer , dimension(:,:,:) :: tgbb
    real(rkx) , pointer , dimension(:,:,:) :: tgrd
    real(rkx) , pointer , dimension(:,:,:) :: tgbrd
    real(rkx) , pointer , dimension(:,:,:) :: taf
    real(rkx) , pointer , dimension(:,:,:) :: tlef
    real(rkx) , pointer , dimension(:,:,:) :: sfice
    real(rkx) , pointer , dimension(:,:,:) :: snag
    real(rkx) , pointer , dimension(:,:,:) :: sncv
    real(rkx) , pointer , dimension(:,:,:) :: scvk
    real(rkx) , pointer , dimension(:,:,:) :: emisv
    real(rkx) , pointer , dimension(:,:,:) :: sent
    real(rkx) , pointer , dimension(:,:,:) :: evpr
    real(rkx) , pointer , dimension(:,:,:) :: deltat
    real(rkx) , pointer , dimension(:,:,:) :: deltaq
    real(rkx) , pointer , dimension(:,:,:) :: drag
    real(rkx) , pointer , dimension(:,:,:) :: prcp
    real(rkx) , pointer , dimension(:,:,:) :: snwm
    real(rkx) , pointer , dimension(:,:,:) :: trnof
    real(rkx) , pointer , dimension(:,:,:) :: sigf
    real(rkx) , pointer , dimension(:,:,:) :: sfcp
    real(rkx) , pointer , dimension(:,:,:) :: srnof
    real(rkx) , pointer , dimension(:,:,:) :: xlai
    real(rkx) , pointer , dimension(:,:,:) :: q2m
    real(rkx) , pointer , dimension(:,:,:) :: t2m
    real(rkx) , pointer , dimension(:,:,:) :: u10m
    real(rkx) , pointer , dimension(:,:,:) :: v10m
    real(rkx) , pointer , dimension(:,:,:) :: taux
    real(rkx) , pointer , dimension(:,:,:) :: tauy
    real(rkx) , pointer , dimension(:,:,:) :: swalb
    real(rkx) , pointer , dimension(:,:,:) :: lwalb
    real(rkx) , pointer , dimension(:,:,:) :: swdiralb
    real(rkx) , pointer , dimension(:,:,:) :: lwdiralb
    real(rkx) , pointer , dimension(:,:,:) :: swdifalb
    real(rkx) , pointer , dimension(:,:,:) :: lwdifalb
    real(rkx) , pointer , dimension(:,:,:) :: wt
    real(rkx) , pointer , dimension(:,:,:) :: eta
    real(rkx) , pointer , dimension(:,:,:) :: hi
    real(rkx) , pointer , dimension(:,:,:) :: hsnow
    real(rkx) , pointer , dimension(:,:,:) :: um10
    real(rkx) , pointer , dimension(:,:,:,:) :: tlake
    logical , pointer , dimension(:,:,:) :: lakmsk
    real(rkx) , pointer , dimension(:,:,:) :: deltas
    real(rkx) , pointer , dimension(:,:,:) :: tdeltas
    real(rkx) , pointer , dimension(:,:,:) :: tskin
    real(rkx) , pointer , dimension(:,:,:) :: sst
    real(rkx) , pointer , dimension(:,:,:) :: zo
    real(rkx) , pointer , dimension(:,:,:) :: ustar
    real(rkx) , pointer , dimension(:,:,:) :: rhoa
#ifdef CLM45
    real(rkx) , pointer , dimension(:,:,:,:) :: vocemiss
    real(rkx) , pointer , dimension(:,:,:,:) :: dustemiss
#endif
  end type lm_state

  type lm_exchange
    real(rkx) , pointer , dimension(:,:) :: ssw2da
    real(rkx) , pointer , dimension(:,:) :: sfracv2d
    real(rkx) , pointer , dimension(:,:) :: sfracb2d
    real(rkx) , pointer , dimension(:,:) :: sfracs2d
    real(rkx) , pointer , dimension(:,:) :: svegfrac2d
    real(rkx) , pointer , dimension(:,:) :: sxlai2d
    real(rkx) , pointer , dimension(:,:,:) :: dailyrnf
    real(rkx) , pointer , dimension(:,:) :: xlat        ! mddom%xlat
    real(rkx) , pointer , dimension(:,:) :: xlon        ! mddom%xlon
    real(rkx) , pointer , dimension(:,:) :: lndcat      ! mddom%lndcat
    real(rkx) , pointer , dimension(:,:) :: ht          ! mddom%ht
    real(rkx) , pointer , dimension(:,:) :: snowam      ! mddom%snowam
    real(rkx) , pointer , dimension(:,:) :: smoist      ! mddom%smoist
    real(rkx) , pointer , dimension(:,:,:) :: rmoist    ! mddom%rmoist
    integer(ik4) , pointer , dimension(:,:) :: iveg     ! mddom%iveg
    integer(ik4) , pointer , dimension(:,:) :: ldmsk    ! mddom%ldmsk
    real(rkx) , pointer , dimension(:,:,:) :: ht1       ! mdsub%ht
    real(rkx) , pointer , dimension(:,:,:) :: lndcat1   ! mdsub%lndcat
    real(rkx) , pointer , dimension(:,:,:) :: xlat1     ! mdsub%xlat
    real(rkx) , pointer , dimension(:,:,:) :: xlon1     ! mdsub%xlon
    real(rkx) , pointer , dimension(:,:,:) :: dhlake1   ! mdsub%dhlake
    integer(ik4) , pointer , dimension(:,:,:) :: ldmsk1 ! mdsub%ldmsk
    integer(ik4) , pointer , dimension(:,:,:) :: iveg1  ! mdsub%iveg
    integer(ik4) , pointer , dimension(:,:) :: icplmsk  ! cplmsk
    real(rkx) , pointer , dimension(:,:) :: patm        ! atms%pb3d(:,:,kz)
    real(rkx) , pointer , dimension(:,:) :: uatm        ! atms%ubx3d(:,:,kz)
    real(rkx) , pointer , dimension(:,:) :: vatm        ! atms%vbx3d(:,:,kz)
    real(rkx) , pointer , dimension(:,:) :: tatm        ! atms%tb3d(:,:,kz)
    real(rkx) , pointer , dimension(:,:) :: thatm       ! atms%th3d(:,:,kz)
    real(rkx) , pointer , dimension(:,:) :: qvatm       ! atms%qxb3d(:,:,kz,iqv)
    real(rkx) , pointer , dimension(:,:) :: hgt         ! za(:,:,kz)
    real(rkx) , pointer , dimension(:,:) :: hpbl        ! zpbl
    real(rkx) , pointer , dimension(:,:) :: hfx         ! sfs%hfx
    real(rkx) , pointer , dimension(:,:) :: qfx         ! sfs%qfx
    real(rkx) , pointer , dimension(:,:) :: tground1    ! sfs%tga
    real(rkx) , pointer , dimension(:,:) :: tground2    ! sfs%tgb
    real(rkx) , pointer , dimension(:,:) :: sfps        ! sfs%psb
    real(rkx) , pointer , dimension(:,:) :: sfta        ! atms%ts2d
    real(rkx) , pointer , dimension(:,:) :: uvdrag      ! sfs%uvdrag
    real(rkx) , pointer , dimension(:,:) :: tgbb        ! sfs%tgbb
    real(rkx) , pointer , dimension(:,:) :: rhox        ! rhox2d
    real(rkx) , pointer , dimension(:,:) :: rswf        ! fsw
    real(rkx) , pointer , dimension(:,:) :: rlwf        ! flw
    real(rkx) , pointer , dimension(:,:) :: dwrlwf      ! flwd
    real(rkx) , pointer , dimension(:,:) :: zencos      ! coszrs
    real(rkx) , pointer , dimension(:,:) :: ncprate     ! pptnc
    real(rkx) , pointer , dimension(:,:) :: cprate      ! cprate
    real(rkx) , pointer , dimension(:,:) :: vegswab     ! sabveg
    real(rkx) , pointer , dimension(:,:) :: lwalb       ! albvl
    real(rkx) , pointer , dimension(:,:) :: swalb       ! albvs
    real(rkx) , pointer , dimension(:,:) :: swdiralb    ! aldirs
    real(rkx) , pointer , dimension(:,:) :: swdifalb    ! aldifs
    real(rkx) , pointer , dimension(:,:) :: lwdiralb    ! aldirl
    real(rkx) , pointer , dimension(:,:) :: lwdifalb    ! aldifl
    real(rkx) , pointer , dimension(:,:) :: swdir       ! solvs
    real(rkx) , pointer , dimension(:,:) :: swdif       ! solvsd
    real(rkx) , pointer , dimension(:,:) :: lwdir       ! solvl
    real(rkx) , pointer , dimension(:,:) :: lwdif       ! solvld
    real(rkx) , pointer , dimension(:,:) :: solinc      ! sinc
    real(rkx) , pointer , dimension(:,:) :: solar       ! solis
    real(rkx) , pointer , dimension(:,:) :: emissivity  ! emiss
    real(rkx) , pointer , dimension(:,:) :: deltaq      ! sdelq
    real(rkx) , pointer , dimension(:,:) :: deltat      ! sdelt
    real(rkx) , pointer , dimension(:,:,:) :: drydepflx   ! drydepflx
    real(rkx) , pointer , dimension(:,:,:) :: wetdepflx   ! wetdepflx
    integer(ik4) , pointer , dimension(:) :: idust        ! dust indices
    real(rkx) , pointer , dimension(:,:) :: zo          ! zo
    real(rkx) , pointer , dimension(:,:) :: ustar       ! ustar
    real(rkx) , pointer , dimension(:,:) :: rhoa        ! xdens
#ifdef CLM
    real(rkx) , pointer , dimension(:,:,:) :: dep_vels
    real(rkx) , pointer , dimension(:,:) :: voc_em0
    real(rkx) , pointer , dimension(:,:) :: voc_em1
    real(rkx) , pointer , dimension(:,:) :: voc_em2
#endif
  end type lm_exchange

  type mod_2_rad
    real(rkx) , pointer , dimension(:,:,:) :: cldfrc
    real(rkx) , pointer , dimension(:,:,:) :: cldlwc
    real(rkx) , pointer , dimension(:,:,:) :: tatms         ! atms%tb3d
    real(rkx) , pointer , dimension(:,:,:) :: rhatms        ! atms%rhb3d
    real(rkx) , pointer , dimension(:,:,:) :: phatms        ! atms%pb3d
    real(rkx) , pointer , dimension(:,:,:) :: pfatms        ! atms%pf3d
    real(rkx) , pointer , dimension(:,:,:) :: deltaz        ! atms%dzq
    real(rkx) , pointer , dimension(:,:) :: psatms          ! atms%ps2d
    real(rkx) , pointer , dimension(:,:,:,:) :: qxatms      ! atms%qxb3d
    real(rkx) , pointer , dimension(:,:,:,:) :: chiatms     ! atms%chib3d
    real(rkx) , pointer , dimension(:,:) :: tg              ! sfs%tgbb
    real(rkx) , pointer , dimension(:,:) :: xlat            ! mddom%xlat
    real(rkx) , pointer , dimension(:,:) :: xlon            ! mddom%xlon
    real(rkx) , pointer , dimension(:,:) :: ptrop
    real(rkx) , pointer , dimension(:,:) :: coszrs
    real(rkx) , pointer , dimension(:,:) :: albvs
    real(rkx) , pointer , dimension(:,:) :: albvl
    real(rkx) , pointer , dimension(:,:) :: aldirs
    real(rkx) , pointer , dimension(:,:) :: aldifs
    real(rkx) , pointer , dimension(:,:) :: aldirl
    real(rkx) , pointer , dimension(:,:) :: aldifl
    real(rkx) , pointer , dimension(:,:) :: emiss
    integer(ik4) , pointer , dimension(:,:) :: ldmsk
  end type mod_2_rad

  type rad_2_mod
    real(rkx) , pointer , dimension(:,:) :: solis
    real(rkx) , pointer , dimension(:,:) :: sinc
    real(rkx) , pointer , dimension(:,:) :: sabveg
    real(rkx) , pointer , dimension(:,:) :: sols
    real(rkx) , pointer , dimension(:,:) :: soll
    real(rkx) , pointer , dimension(:,:) :: solvs
    real(rkx) , pointer , dimension(:,:) :: solvsd
    real(rkx) , pointer , dimension(:,:) :: solvl
    real(rkx) , pointer , dimension(:,:) :: solvld
    real(rkx) , pointer , dimension(:,:) :: fsw
    real(rkx) , pointer , dimension(:,:) :: flw
    real(rkx) , pointer , dimension(:,:) :: flwd
    real(rkx) , pointer , dimension(:,:,:) :: heatrt
  end type rad_2_mod

  type mod_2_cum
    real(rkx) , pointer , dimension(:,:) :: ht        ! mddom%ht
    real(rkx) , pointer , dimension(:,:) :: psa       ! sfs%psa
    real(rkx) , pointer , dimension(:,:) :: psb       ! sfs%psb
    real(rkx) , pointer , dimension(:,:) :: psdotb    ! sfs%psdotb
    real(rkx) , pointer , dimension(:,:) :: psf       ! atms%ps2d
    real(rkx) , pointer , dimension(:,:,:) :: pas     ! atms%pb3d
    real(rkx) , pointer , dimension(:,:,:) :: pasf    ! atms%pf3d
    real(rkx) , pointer , dimension(:,:,:) :: zas     ! atms%za
    real(rkx) , pointer , dimension(:,:,:) :: tas     ! atms%tb3d
    real(rkx) , pointer , dimension(:,:,:) :: uas     ! atms%ubx3d
    real(rkx) , pointer , dimension(:,:,:) :: vas     ! atms%vbx3d
    real(rkx) , pointer , dimension(:,:,:) :: was     ! atms%wx3d
    real(rkx) , pointer , dimension(:,:,:) :: wpas    ! atms%wpx3d
    real(rkx) , pointer , dimension(:,:,:) :: qsas    ! atms%qsb3d
    real(rkx) , pointer , dimension(:,:,:) :: tkeas   ! atms%tke
    real(rkx) , pointer , dimension(:,:,:) :: rhoas   ! atms%rhob3d
    real(rkx) , pointer , dimension(:,:,:) :: zfs     ! atms%zq
    real(rkx) , pointer , dimension(:,:,:) :: dzq     ! atms%dzq
    real(rkx) , pointer , dimension(:,:,:) :: qq1     ! atm1%q
    real(rkx) , pointer , dimension(:,:,:) :: qdot    ! qdot
    real(rkx) , pointer , dimension(:,:,:,:) :: qxas  ! atms%qxb3d
    real(rkx) , pointer , dimension(:,:,:,:) :: chias ! atms%chib3d
    real(rkx) , pointer , dimension(:,:) :: qfx       ! sfs%qfx
    real(rkx) , pointer , dimension(:,:) :: hfx       ! sfs%hfx
    real(rkx) , pointer , dimension(:,:,:) :: ccn     ! ccn
    integer(ik4) , pointer , dimension(:,:) :: ktrop
    integer(ik4) , pointer , dimension(:,:) :: ldmsk
  end type mod_2_cum

  type cum_2_mod
    real(rkx) , pointer , dimension(:,:,:) :: tten     ! aten%t
    real(rkx) , pointer , dimension(:,:,:) :: uten     ! aten%u
    real(rkx) , pointer , dimension(:,:,:) :: vten     ! aten%v
    real(rkx) , pointer , dimension(:,:,:,:) :: qxten  ! aten%qx
    real(rkx) , pointer , dimension(:,:,:,:) :: chiten ! chiten
    real(rkx) , pointer , dimension(:,:) :: rainc
    real(rkx) , pointer , dimension(:,:) :: pcratec
    real(rkx) , pointer , dimension(:,:,:) :: convpr
    real(rkx) , pointer , dimension(:,:,:) :: cldfrc
    real(rkx) , pointer , dimension(:,:,:) :: cldlwc
    real(rkx) , pointer , dimension(:,:,:) :: q_detr
    real(rkx) , pointer , dimension(:,:,:) :: rain_cc
    integer(ik4) , pointer , dimension(:,:) :: kcumtop
    integer(ik4) , pointer , dimension(:,:) :: kcumbot
  end type cum_2_mod

  type mod_2_pbl
    real(rkx) , pointer , dimension(:,:) :: coriol      ! mddom%coriol
    real(rkx) , pointer , dimension(:,:) :: psdot       ! psdot
    real(rkx) , pointer , dimension(:,:) :: psb         ! sfs%psb
    real(rkx) , pointer , dimension(:,:) :: tgb         ! sfs%tgb
    real(rkx) , pointer , dimension(:,:) :: qfx         ! sfs%qfx
    real(rkx) , pointer , dimension(:,:) :: hfx         ! sfs%hfx
    real(rkx) , pointer , dimension(:,:) :: uvdrag      ! sfs%uvdrag
    real(rkx) , pointer , dimension(:,:,:) :: uxatm     ! atms%ubx3d
    real(rkx) , pointer , dimension(:,:,:) :: vxatm     ! atms%vbx3d
    real(rkx) , pointer , dimension(:,:,:) :: udatm     ! atms%ubd3d
    real(rkx) , pointer , dimension(:,:,:) :: vdatm     ! atms%vbd3d
    real(rkx) , pointer , dimension(:,:,:) :: tatm      ! atms%tb3d
    real(rkx) , pointer , dimension(:,:,:) :: patm      ! atms%pb3d
    real(rkx) , pointer , dimension(:,:,:) :: patmf     ! atms%pf3d
    real(rkx) , pointer , dimension(:,:,:,:) :: qxatm   ! atms%qx
    real(rkx) , pointer , dimension(:,:,:) :: tkests    ! atms%tke
    real(rkx) , pointer , dimension(:,:,:) :: thatm     ! atms%th3d
    real(rkx) , pointer , dimension(:,:,:) :: tpatm     ! atms%tp3d
    real(rkx) , pointer , dimension(:,:,:) :: za        ! atms%za
    real(rkx) , pointer , dimension(:,:,:) :: zq        ! atms%zq
    real(rkx) , pointer , dimension(:,:,:) :: dzq       ! atms%dzq
    real(rkx) , pointer , dimension(:,:) :: rhox2d      ! atms%rhox2d
    real(rkx) , pointer , dimension(:,:,:) :: heatrt    ! heatrt
    real(rkx) , pointer , dimension(:,:,:,:) :: chib    ! chib
    real(rkx) , pointer , dimension(:,:,:) :: chifxuw   ! chifxuw
    real(rkx) , pointer , dimension(:,:,:) :: drydepv   ! drydepv
    integer(ik4) , pointer , dimension(:,:) :: ktrop    ! ktrop
  end type mod_2_pbl

  type pbl_2_mod
    real(rkx) , pointer , dimension(:,:,:) :: tten       ! aten%t
    real(rkx) , pointer , dimension(:,:,:) :: uten       ! aten%u
    real(rkx) , pointer , dimension(:,:,:) :: vten       ! aten%v
    real(rkx) , pointer , dimension(:,:,:,:) :: qxten    ! aten%qx
    real(rkx) , pointer , dimension(:,:,:) :: tketen     ! aten%tke
    real(rkx) , pointer , dimension(:,:,:) :: uuwten     ! uwten%u
    real(rkx) , pointer , dimension(:,:,:) :: vuwten     ! uwten%v
    real(rkx) , pointer , dimension(:,:,:) :: tuwten     ! uwten%t
    real(rkx) , pointer , dimension(:,:,:) :: tkeuwten   ! uwten%tke
    real(rkx) , pointer , dimension(:,:,:,:) :: qxuwten  ! uwten%qx
    real(rkx) , pointer , dimension(:,:,:) :: difft      ! adf%difft
    real(rkx) , pointer , dimension(:,:,:,:) :: diffqx   ! adf%diffqx
    real(rkx) , pointer , dimension(:,:,:,:) :: diagqx   ! holtten%qx
    real(rkx) , pointer , dimension(:,:,:,:) :: chiten   ! chiten
    real(rkx) , pointer , dimension(:,:,:) :: remdrd     ! remdrd
    real(rkx) , pointer , dimension(:,:) :: zpbl
    integer(ik4) , pointer , dimension(:,:) :: kpbl
  end type pbl_2_mod

end module mod_regcm_types

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
