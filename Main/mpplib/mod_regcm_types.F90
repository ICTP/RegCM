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
    integer(ik4) :: ibt1 , ibt2 , ibt4 , ibb1 , ibb2 , ibb4
    integer(ik4) :: jbl1 , jbl2 , jbl4 , jbr1 , jbr2 , jbr4
  end type model_area

  type domain
    real(rk8) , pointer , dimension(:,:) :: ht
    real(rk8) , pointer , dimension(:,:) :: lndcat
    real(rk8) , pointer , dimension(:,:) :: xlat
    real(rk8) , pointer , dimension(:,:) :: xlon
    real(rk8) , pointer , dimension(:,:) :: mask
    real(rk8) , pointer , dimension(:,:) :: dlat
    real(rk8) , pointer , dimension(:,:) :: dlon
    real(rk8) , pointer , dimension(:,:) :: msfx
    real(rk8) , pointer , dimension(:,:) :: msfd
    real(rk8) , pointer , dimension(:,:) :: coriol
    real(rk8) , pointer , dimension(:,:) :: ef
    real(rk8) , pointer , dimension(:,:) :: ddx
    real(rk8) , pointer , dimension(:,:) :: ddy
    real(rk8) , pointer , dimension(:,:) :: ex
    real(rk8) , pointer , dimension(:,:) :: crx
    real(rk8) , pointer , dimension(:,:) :: cry
    real(rk8) , pointer , dimension(:,:) :: dmdy
    real(rk8) , pointer , dimension(:,:) :: dmdx
    real(rk8) , pointer , dimension(:,:) :: snowam
    real(rk8) , pointer , dimension(:,:) :: satmoist
    real(rk8) , pointer , dimension(:,:) :: dhlake
    integer(ik4) , pointer , dimension(:,:) :: ldmsk
    integer(ik4) , pointer , dimension(:,:) :: iveg
  end type domain

  type domain_subgrid
    real(rk8) , pointer , dimension(:,:,:) :: ht
    real(rk8) , pointer , dimension(:,:,:) :: lndcat
    real(rk8) , pointer , dimension(:,:,:) :: xlat
    real(rk8) , pointer , dimension(:,:,:) :: xlon
    real(rk8) , pointer , dimension(:,:,:) :: mask
    real(rk8) , pointer , dimension(:,:,:) :: dhlake
    integer(ik4) , pointer , dimension(:,:,:) :: ldmsk
    integer(ik4) , pointer , dimension(:,:,:) :: iveg
  end type domain_subgrid

  type reference_atmosphere
    real(rk8) , pointer , dimension(:,:) :: ps
    real(rk8) , pointer , dimension(:,:,:) :: t
    real(rk8) , pointer , dimension(:,:,:) :: pr
    real(rk8) , pointer , dimension(:,:,:) :: rho
  end type reference_atmosphere

  type atmstate_a
    real(rk8) , pointer , dimension(:,:,:) :: u
    real(rk8) , pointer , dimension(:,:,:) :: v
    real(rk8) , pointer , dimension(:,:,:) :: w
    real(rk8) , pointer , dimension(:,:,:) :: t
    real(rk8) , pointer , dimension(:,:,:,:) :: qx
    real(rk8) , pointer , dimension(:,:,:) :: tke
    real(rk8) , pointer , dimension(:,:,:) :: pp
    real(rk8) , pointer , dimension(:,:,:) :: pr
    real(rk8) , pointer , dimension(:,:,:) :: rho
  end type atmstate_a

  type atmstate_b
    real(rk8) , pointer , dimension(:,:,:) :: u
    real(rk8) , pointer , dimension(:,:,:) :: v
    real(rk8) , pointer , dimension(:,:,:) :: w
    real(rk8) , pointer , dimension(:,:,:) :: t
    real(rk8) , pointer , dimension(:,:,:,:) :: qx
    real(rk8) , pointer , dimension(:,:,:) :: tke
    real(rk8) , pointer , dimension(:,:,:) :: pp
    real(rk8) , pointer , dimension(:,:,:) :: pr
    real(rk8) , pointer , dimension(:,:,:) :: rho
  end type atmstate_b

  type atmstate_c
    real(rk8) , pointer , dimension(:,:,:) :: u
    real(rk8) , pointer , dimension(:,:,:) :: v
    real(rk8) , pointer , dimension(:,:,:) :: w
    real(rk8) , pointer , dimension(:,:,:) :: t
    real(rk8) , pointer , dimension(:,:,:,:) :: qx
    real(rk8) , pointer , dimension(:,:,:) :: tke
    real(rk8) , pointer , dimension(:,:,:) :: pp
    real(rk8) , pointer , dimension(:,:,:) :: pr
    real(rk8) , pointer , dimension(:,:,:) :: rho
  end type atmstate_c

  type atmstate_tendency
    real(rk8) , pointer , dimension(:,:,:) :: u
    real(rk8) , pointer , dimension(:,:,:) :: v
    real(rk8) , pointer , dimension(:,:,:) :: w
    real(rk8) , pointer , dimension(:,:,:) :: t
    real(rk8) , pointer , dimension(:,:,:,:) :: qx
    real(rk8) , pointer , dimension(:,:,:) :: tke
    real(rk8) , pointer , dimension(:,:,:) :: pp
    real(rk8) , pointer , dimension(:,:,:) :: pr
    real(rk8) , pointer , dimension(:,:,:) :: rho
  end type atmstate_tendency

  type atmstate_decoupled
    real(rk8) , pointer , dimension(:,:,:) :: u
    real(rk8) , pointer , dimension(:,:,:) :: v
    real(rk8) , pointer , dimension(:,:,:) :: w
    real(rk8) , pointer , dimension(:,:,:) :: t
    real(rk8) , pointer , dimension(:,:,:,:) :: qx
    real(rk8) , pointer , dimension(:,:,:) :: tke
    real(rk8) , pointer , dimension(:,:,:) :: pp
    real(rk8) , pointer , dimension(:,:,:) :: pr
    real(rk8) , pointer , dimension(:,:,:) :: rho
  end type atmstate_decoupled

  type tcm_state
    ! TKE*ps
    real(rk8) , pointer , dimension(:,:,:) :: tkeps  ! (m^2/s^2 * cb)
    ! Coupled TKE Advective Tendency
    real(rk8) , pointer , dimension(:,:,:) :: advtke ! (m^2/s^3 * cb)
    ! Vertical momentum diffusivity
    real(rk8) , pointer , dimension(:,:,:) :: kzm    ! (m^2/s)
    ! Vertical scalar diffusivity
    real(rk8) , pointer , dimension(:,:,:) :: kth    ! (m^2/s)
    ! Boundary layer height (m)
    real(rk8) , pointer , dimension(:,:) :: zpbl     ! (m)
    ! Surface layer TKE
    real(rk8) , pointer , dimension(:,:) :: srftke   ! (m^2/s^2)
  end type tcm_state

  type tendiag
    real(rk8) , pointer , dimension(:,:,:) :: adh
    real(rk8) , pointer , dimension(:,:,:) :: adv
    real(rk8) , pointer , dimension(:,:,:) :: tbl
    real(rk8) , pointer , dimension(:,:,:) :: dif
    real(rk8) , pointer , dimension(:,:,:) :: bdy
    real(rk8) , pointer , dimension(:,:,:) :: con
    real(rk8) , pointer , dimension(:,:,:) :: adi
    real(rk8) , pointer , dimension(:,:,:) :: rad
    real(rk8) , pointer , dimension(:,:,:) :: lsc
  end type tendiag

  type qendiag
    real(rk8) , pointer , dimension(:,:,:) :: adh
    real(rk8) , pointer , dimension(:,:,:) :: adv
    real(rk8) , pointer , dimension(:,:,:) :: tbl
    real(rk8) , pointer , dimension(:,:,:) :: dif
    real(rk8) , pointer , dimension(:,:,:) :: bdy
    real(rk8) , pointer , dimension(:,:,:) :: con
    real(rk8) , pointer , dimension(:,:,:) :: adi
    real(rk8) , pointer , dimension(:,:,:) :: rad
    real(rk8) , pointer , dimension(:,:,:) :: lsc
  end type qendiag

  type surfstate
    real(rk8) , pointer , dimension(:,:) :: psa
    real(rk8) , pointer , dimension(:,:) :: psb
    real(rk8) , pointer , dimension(:,:) :: psc
    real(rk8) , pointer , dimension(:,:) :: psdota
    real(rk8) , pointer , dimension(:,:) :: psdotb
    real(rk8) , pointer , dimension(:,:) :: tga
    real(rk8) , pointer , dimension(:,:) :: tgb
    real(rk8) , pointer , dimension(:,:) :: rainc
    real(rk8) , pointer , dimension(:,:) :: rainnc
    real(rk8) , pointer , dimension(:,:) :: snownc
    real(rk8) , pointer , dimension(:,:) :: hfx
    real(rk8) , pointer , dimension(:,:) :: qfx
    real(rk8) , pointer , dimension(:,:) :: tgbb
    real(rk8) , pointer , dimension(:,:) :: uvdrag
  end type surfstate

  type slice
    real(rk8) , pointer , dimension(:,:,:) :: tb3d
    real(rk8) , pointer , dimension(:,:,:) :: thx3d
    real(rk8) , pointer , dimension(:,:,:) :: pb3d
    real(rk8) , pointer , dimension(:,:,:) :: pf3d
    real(rk8) , pointer , dimension(:,:,:) :: rhob3d
    real(rk8) , pointer , dimension(:,:,:) :: ubx3d
    real(rk8) , pointer , dimension(:,:,:) :: vbx3d
    real(rk8) , pointer , dimension(:,:,:) :: wpx3d
    real(rk8) , pointer , dimension(:,:,:) :: ubd3d
    real(rk8) , pointer , dimension(:,:,:) :: vbd3d
    real(rk8) , pointer , dimension(:,:,:) :: rhb3d
    real(rk8) , pointer , dimension(:,:,:) :: qsb3d
    real(rk8) , pointer , dimension(:,:,:,:) :: qxb3d
    real(rk8) , pointer , dimension(:,:,:) :: zq
    real(rk8) , pointer , dimension(:,:,:) :: za
    real(rk8) , pointer , dimension(:,:,:) :: dzq
    real(rk8) , pointer , dimension(:,:) :: rhox2d
    real(rk8) , pointer , dimension(:,:) :: ps2d
    real(rk8) , pointer , dimension(:,:,:) :: tkeb3d
    real(rk8) , pointer , dimension(:,:,:,:) :: chib3d
  end type slice

  type diffx
    real(rk8) , pointer , dimension(:,:,:) :: difft
    real(rk8) , pointer , dimension(:,:,:) :: difuu
    real(rk8) , pointer , dimension(:,:,:) :: difuv
    real(rk8) , pointer , dimension(:,:,:) :: diffw
    real(rk8) , pointer , dimension(:,:,:,:) :: diffqx
  end type diffx

  type v3dbound
    real(rk8) , pointer , dimension(:,:,:) :: b0
    real(rk8) , pointer , dimension(:,:,:) :: b1
    real(rk8) , pointer , dimension(:,:,:) :: bt
  end type v3dbound

  type v2dbound
    real(rk8) , pointer , dimension(:,:) :: b0
    real(rk8) , pointer , dimension(:,:) :: b1
    real(rk8) , pointer , dimension(:,:) :: bt
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
    real(rk8) , pointer , dimension(:,:) :: psfc
    real(rk8) , pointer , dimension(:,:) :: tsfc
    real(rk8) , pointer , dimension(:,:) :: qsfc
    real(rk8) , pointer , dimension(:,:) :: swrd
    real(rk8) , pointer , dimension(:,:) :: lwrd
    real(rk8) , pointer , dimension(:,:) :: dlwr
    real(rk8) , pointer , dimension(:,:) :: lhfx
    real(rk8) , pointer , dimension(:,:) :: shfx
    real(rk8) , pointer , dimension(:,:) :: prec
    real(rk8) , pointer , dimension(:,:) :: wndu
    real(rk8) , pointer , dimension(:,:) :: wndv
    real(rk8) , pointer , dimension(:,:) :: rnof
    real(rk8) , pointer , dimension(:,:) :: snof
    real(rk8) , pointer , dimension(:,:) :: taux
    real(rk8) , pointer , dimension(:,:) :: tauy
    real(rk8) , pointer , dimension(:,:) :: wspd
    real(rk8) , pointer , dimension(:,:) :: nflx
    real(rk8) , pointer , dimension(:,:) :: sflx
    real(rk8) , pointer , dimension(:,:) :: snow
    real(rk8) , pointer , dimension(:,:) :: dswr
  end type exp_data

  type imp_data
    real(rk8) , pointer , dimension(:,:) :: sst
    real(rk8) , pointer , dimension(:,:) :: sit
    real(rk8) , pointer , dimension(:,:) :: msk
  end type imp_data

  type lm_state
    real(rk8) , pointer , dimension(:,:,:) :: gwet
    real(rk8) , pointer , dimension(:,:,:) :: ssw
    real(rk8) , pointer , dimension(:,:,:) :: rsw
    real(rk8) , pointer , dimension(:,:,:) :: tsw
    real(rk8) , pointer , dimension(:,:,:) :: ldew
    real(rk8) , pointer , dimension(:,:,:) :: lncl
    real(rk8) , pointer , dimension(:,:,:) :: tgbb
    real(rk8) , pointer , dimension(:,:,:) :: tgrd
    real(rk8) , pointer , dimension(:,:,:) :: tgbrd
    real(rk8) , pointer , dimension(:,:,:) :: taf
    real(rk8) , pointer , dimension(:,:,:) :: tlef
    real(rk8) , pointer , dimension(:,:,:) :: sfice
    real(rk8) , pointer , dimension(:,:,:) :: snag
    real(rk8) , pointer , dimension(:,:,:) :: sncv
    real(rk8) , pointer , dimension(:,:,:) :: scvk
    real(rk8) , pointer , dimension(:,:,:) :: emisv
    real(rk8) , pointer , dimension(:,:,:) :: sent
    real(rk8) , pointer , dimension(:,:,:) :: evpr
    real(rk8) , pointer , dimension(:,:,:) :: deltat
    real(rk8) , pointer , dimension(:,:,:) :: deltaq
    real(rk8) , pointer , dimension(:,:,:) :: drag
    real(rk8) , pointer , dimension(:,:,:) :: prcp
    real(rk8) , pointer , dimension(:,:,:) :: snwm
    real(rk8) , pointer , dimension(:,:,:) :: trnof
    real(rk8) , pointer , dimension(:,:,:) :: sigf
    real(rk8) , pointer , dimension(:,:,:) :: sfcp
    real(rk8) , pointer , dimension(:,:,:) :: srnof
    real(rk8) , pointer , dimension(:,:,:) :: xlai
    real(rk8) , pointer , dimension(:,:,:) :: q2m
    real(rk8) , pointer , dimension(:,:,:) :: t2m
    real(rk8) , pointer , dimension(:,:,:) :: u10m
    real(rk8) , pointer , dimension(:,:,:) :: v10m
    real(rk8) , pointer , dimension(:,:,:) :: taux
    real(rk8) , pointer , dimension(:,:,:) :: tauy
    real(rk8) , pointer , dimension(:,:,:) :: swalb
    real(rk8) , pointer , dimension(:,:,:) :: lwalb
    real(rk8) , pointer , dimension(:,:,:) :: swdiralb
    real(rk8) , pointer , dimension(:,:,:) :: lwdiralb
    real(rk8) , pointer , dimension(:,:,:) :: swdifalb
    real(rk8) , pointer , dimension(:,:,:) :: lwdifalb
    real(rk8) , pointer , dimension(:,:,:) :: wt
    real(rk8) , pointer , dimension(:,:,:) :: eta
    real(rk8) , pointer , dimension(:,:,:) :: hi
    real(rk8) , pointer , dimension(:,:,:) :: hsnow
    real(rk8) , pointer , dimension(:,:,:,:) :: tlake
    logical , pointer , dimension(:,:,:) :: lakmsk
    real(rk8) , pointer , dimension(:,:,:) :: deltas
    real(rk8) , pointer , dimension(:,:,:) :: tdeltas
    real(rk8) , pointer , dimension(:,:,:) :: tskin
    real(rk8) , pointer , dimension(:,:,:) :: sst
    real(rk8) , pointer , dimension(:,:,:,:) :: vocemiss
  end type lm_state

  type lm_exchange
    real(rk8) , pointer , dimension(:,:) :: ssw2da
    real(rk8) , pointer , dimension(:,:) :: sfracv2d
    real(rk8) , pointer , dimension(:,:) :: sfracb2d
    real(rk8) , pointer , dimension(:,:) :: sfracs2d
    real(rk8) , pointer , dimension(:,:) :: svegfrac2d
    real(rk8) , pointer , dimension(:,:) :: sxlai2d
    real(rk8) , pointer , dimension(:,:,:) :: dailyrnf
    real(rk8) , pointer , dimension(:,:) :: xlat        ! mddom%xlat
    real(rk8) , pointer , dimension(:,:) :: xlon        ! mddom%xlon
    real(rk8) , pointer , dimension(:,:) :: lndcat      ! mddom%lndcat
    real(rk8) , pointer , dimension(:,:) :: ht          ! mddom%ht
    real(rk8) , pointer , dimension(:,:) :: snowam      ! mddom%snowam
    real(rk8) , pointer , dimension(:,:) :: smoist      ! mddom%satmoist
    integer(ik4) , pointer , dimension(:,:) :: iveg     ! mddom%iveg
    integer(ik4) , pointer , dimension(:,:) :: ldmsk    ! mddom%ldmsk
    real(rk8) , pointer , dimension(:,:,:) :: ht1       ! mdsub%ht
    real(rk8) , pointer , dimension(:,:,:) :: lndcat1   ! mdsub%lndcat
    real(rk8) , pointer , dimension(:,:,:) :: xlat1     ! mdsub%xlat
    real(rk8) , pointer , dimension(:,:,:) :: xlon1     ! mdsub%xlon
    real(rk8) , pointer , dimension(:,:,:) :: dhlake1   ! mdsub%dhlake
    integer(ik4) , pointer , dimension(:,:,:) :: ldmsk1 ! mdsub%ldmsk
    integer(ik4) , pointer , dimension(:,:,:) :: iveg1  ! mdsub%iveg
    integer(ik4) , pointer , dimension(:,:) :: icplmsk  ! cplmsk
    real(rk8) , pointer , dimension(:,:) :: patm        ! atms%pb3d(:,:,kz)
    real(rk8) , pointer , dimension(:,:) :: uatm        ! atms%ubx3d(:,:,kz)
    real(rk8) , pointer , dimension(:,:) :: vatm        ! atms%vbx3d(:,:,kz)
    real(rk8) , pointer , dimension(:,:) :: tatm        ! atms%tb3d(:,:,kz)
    real(rk8) , pointer , dimension(:,:) :: thatm       ! atms%thx3d(:,:,kz)
    real(rk8) , pointer , dimension(:,:) :: qvatm       ! atms%qxb3d(:,:,kz,iqv)
    real(rk8) , pointer , dimension(:,:) :: hgt         ! za(:,:,kz)
    real(rk8) , pointer , dimension(:,:) :: hpbl        ! zpbl
    real(rk8) , pointer , dimension(:,:) :: hfx         ! sfs%hfx
    real(rk8) , pointer , dimension(:,:) :: qfx         ! sfs%qfx
    real(rk8) , pointer , dimension(:,:) :: tground1    ! sfs%tga
    real(rk8) , pointer , dimension(:,:) :: tground2    ! sfs%tgb
    real(rk8) , pointer , dimension(:,:) :: sfps        ! sfs%psb
    real(rk8) , pointer , dimension(:,:) :: uvdrag      ! sfs%uvdrag
    real(rk8) , pointer , dimension(:,:) :: tgbb        ! sfs%tgbb
    real(rk8) , pointer , dimension(:,:) :: rhox        ! rhox2d
    real(rk8) , pointer , dimension(:,:) :: rswf        ! fsw
    real(rk8) , pointer , dimension(:,:) :: rlwf        ! flw
    real(rk8) , pointer , dimension(:,:) :: dwrlwf      ! flwd
    real(rk8) , pointer , dimension(:,:) :: zencos      ! coszrs
    real(rk8) , pointer , dimension(:,:) :: ncprate     ! pptnc
    real(rk8) , pointer , dimension(:,:) :: cprate      ! cprate
    real(rk8) , pointer , dimension(:,:) :: vegswab     ! sabveg
    real(rk8) , pointer , dimension(:,:) :: lwalb       ! albvl
    real(rk8) , pointer , dimension(:,:) :: swalb       ! albvs
    real(rk8) , pointer , dimension(:,:) :: swdiralb    ! aldirs
    real(rk8) , pointer , dimension(:,:) :: swdifalb    ! aldifs
    real(rk8) , pointer , dimension(:,:) :: lwdiralb    ! aldirl
    real(rk8) , pointer , dimension(:,:) :: lwdifalb    ! aldifl
    real(rk8) , pointer , dimension(:,:) :: swdir       ! solvs
    real(rk8) , pointer , dimension(:,:) :: swdif       ! solvsd
    real(rk8) , pointer , dimension(:,:) :: lwdir       ! solvl
    real(rk8) , pointer , dimension(:,:) :: lwdif       ! solvld
    real(rk8) , pointer , dimension(:,:) :: solinc      ! sinc
    real(rk8) , pointer , dimension(:,:) :: solar       ! solis
    real(rk8) , pointer , dimension(:,:) :: emissivity  ! emiss
    real(rk8) , pointer , dimension(:,:) :: deltaq      ! sdelq
    real(rk8) , pointer , dimension(:,:) :: deltat      ! sdelt
  end type lm_exchange

  type mod_2_rad
    real(rk8) , pointer , dimension(:,:) :: psb             ! sfs%psb
    real(rk8) , pointer , dimension(:,:,:) :: cldfrc
    real(rk8) , pointer , dimension(:,:,:) :: cldlwc
    real(rk8) , pointer , dimension(:,:,:) :: tatms         ! atms%tb3d
    real(rk8) , pointer , dimension(:,:,:) :: rhatms        ! atms%rhb3d
    real(rk8) , pointer , dimension(:,:,:) :: phatms        ! atms%pb3d
    real(rk8) , pointer , dimension(:,:,:) :: pfatms        ! atms%pf3d
    real(rk8) , pointer , dimension(:,:) :: psatms          ! atms%ps2d
    real(rk8) , pointer , dimension(:,:,:,:) :: qxatms      ! atms%qxb3d
    real(rk8) , pointer , dimension(:,:,:,:) :: chiatms     ! atms%chib3d
    real(rk8) , pointer , dimension(:,:) :: tg              ! sfs%tgbb
    real(rk8) , pointer , dimension(:,:) :: xlat            ! mddom%xlat
    real(rk8) , pointer , dimension(:,:) :: xlon            ! mddom%xlon
    real(rk8) , pointer , dimension(:,:) :: ptrop
    real(rk8) , pointer , dimension(:,:) :: coszrs
    real(rk8) , pointer , dimension(:,:) :: albvs
    real(rk8) , pointer , dimension(:,:) :: albvl
    real(rk8) , pointer , dimension(:,:) :: aldirs
    real(rk8) , pointer , dimension(:,:) :: aldifs
    real(rk8) , pointer , dimension(:,:) :: aldirl
    real(rk8) , pointer , dimension(:,:) :: aldifl
    real(rk8) , pointer , dimension(:,:) :: emiss
    integer(ik4) , pointer , dimension(:,:) :: ldmsk
  end type mod_2_rad

  type rad_2_mod
    real(rk8) , pointer , dimension(:,:) :: solis
    real(rk8) , pointer , dimension(:,:) :: sinc
    real(rk8) , pointer , dimension(:,:) :: sabveg
    real(rk8) , pointer , dimension(:,:) :: sols
    real(rk8) , pointer , dimension(:,:) :: soll
    real(rk8) , pointer , dimension(:,:) :: solvs
    real(rk8) , pointer , dimension(:,:) :: solvsd
    real(rk8) , pointer , dimension(:,:) :: solvl
    real(rk8) , pointer , dimension(:,:) :: solvld
    real(rk8) , pointer , dimension(:,:) :: fsw
    real(rk8) , pointer , dimension(:,:) :: flw
    real(rk8) , pointer , dimension(:,:) :: flwd
    real(rk8) , pointer , dimension(:,:,:) :: heatrt
  end type rad_2_mod

  type mod_2_cum
    real(rk8) , pointer , dimension(:,:) :: ht        ! mddom%ht
    real(rk8) , pointer , dimension(:,:) :: psb       ! sfs%psb
    real(rk8) , pointer , dimension(:,:,:) :: pas     ! atms%pb3d
    real(rk8) , pointer , dimension(:,:,:) :: pasf    ! atms%pf3d
    real(rk8) , pointer , dimension(:,:,:) :: zas     ! atms%za
    real(rk8) , pointer , dimension(:,:,:) :: tas     ! atms%tb3d
    real(rk8) , pointer , dimension(:,:,:) :: uas     ! atms%ubx3d
    real(rk8) , pointer , dimension(:,:,:) :: vas     ! atms%vbx3d
    real(rk8) , pointer , dimension(:,:,:) :: wpas    ! atms%wpx3d
    real(rk8) , pointer , dimension(:,:,:) :: qsas    ! atms%qsb3d
    real(rk8) , pointer , dimension(:,:,:) :: rhoas   ! atms%rhob3d
    real(rk8) , pointer , dimension(:,:,:) :: zfs     ! atms%zq
    real(rk8) , pointer , dimension(:,:,:) :: dzq     ! atms%dzq
    real(rk8) , pointer , dimension(:,:,:) :: qdot    ! qdot
    real(rk8) , pointer , dimension(:,:,:,:) :: qxas  ! atms%qxb3d
    real(rk8) , pointer , dimension(:,:,:,:) :: chias ! atms%chib3d
    real(rk8) , pointer , dimension(:,:) :: qfx       ! sfs%qfx
    real(rk8) , pointer , dimension(:,:) :: hfx       ! sfs%hfx
    integer(ik4) , pointer , dimension(:,:) :: ktrop
    integer(ik4) , pointer , dimension(:,:) :: ldmsk
  end type mod_2_cum

  type cum_2_mod
    real(rk8) , pointer , dimension(:,:,:) :: tten     ! aten%t
    real(rk8) , pointer , dimension(:,:,:) :: uten     ! aten%u
    real(rk8) , pointer , dimension(:,:,:) :: vten     ! aten%v
    real(rk8) , pointer , dimension(:,:,:,:) :: qxten  ! aten%qx
    real(rk8) , pointer , dimension(:,:,:,:) :: chiten ! chiten
    real(rk8) , pointer , dimension(:,:) :: rainc
    real(rk8) , pointer , dimension(:,:) :: pcratec
    real(rk8) , pointer , dimension(:,:,:) :: convpr
    real(rk8) , pointer , dimension(:,:,:) :: cldfrc
    real(rk8) , pointer , dimension(:,:,:) :: cldlwc
    real(rk8) , pointer , dimension(:,:,:) :: q_detr
    integer(ik4) , pointer , dimension(:,:) :: kcumtop
    integer(ik4) , pointer , dimension(:,:) :: kcumbot
  end type cum_2_mod

  type mod_2_pbl
    real(rk8) , pointer , dimension(:,:) :: coriol      ! mddom%coriol
    real(rk8) , pointer , dimension(:,:) :: psdot       ! psdot
    real(rk8) , pointer , dimension(:,:) :: psb         ! sfs%psb
    real(rk8) , pointer , dimension(:,:) :: tgb         ! sfs%tgb
    real(rk8) , pointer , dimension(:,:) :: qfx         ! sfs%qfx
    real(rk8) , pointer , dimension(:,:) :: hfx         ! sfs%hfx
    real(rk8) , pointer , dimension(:,:) :: uvdrag      ! sfs%uvdrag
    real(rk8) , pointer , dimension(:,:,:) :: uxatm     ! atms%ubx3d
    real(rk8) , pointer , dimension(:,:,:) :: vxatm     ! atms%vbx3d
    real(rk8) , pointer , dimension(:,:,:) :: udatm     ! atms%ubd3d
    real(rk8) , pointer , dimension(:,:,:) :: vdatm     ! atms%vbd3d
    real(rk8) , pointer , dimension(:,:,:) :: tatm      ! atms%tb3d
    real(rk8) , pointer , dimension(:,:,:) :: patm      ! atms%pb3d
    real(rk8) , pointer , dimension(:,:,:) :: patmf     ! atms%pf3d
    real(rk8) , pointer , dimension(:,:,:,:) :: qxatm   ! atms%qx
    real(rk8) , pointer , dimension(:,:,:) :: tkests    ! atms%tke
    real(rk8) , pointer , dimension(:,:,:) :: thxatm    ! atms%thx3d
    real(rk8) , pointer , dimension(:,:,:) :: za        ! atms%za
    real(rk8) , pointer , dimension(:,:,:) :: zq        ! atms%zq
    real(rk8) , pointer , dimension(:,:,:) :: dzq       ! atms%dzq
    real(rk8) , pointer , dimension(:,:,:) :: heatrt    ! heatrt
    real(rk8) , pointer , dimension(:,:,:,:) :: chib    ! chib
    real(rk8) , pointer , dimension(:,:,:) :: chifxuw   ! chifxuw
    real(rk8) , pointer , dimension(:,:,:) :: drydepv   ! drydepv
    real(rk8) , pointer , dimension(:,:) :: rhox2d      ! rhox2d
    integer(ik4) , pointer , dimension(:,:) :: ktrop    ! ktrop
  end type mod_2_pbl

  type pbl_2_mod
    real(rk8) , pointer , dimension(:,:,:) :: tten       ! aten%t
    real(rk8) , pointer , dimension(:,:,:) :: uten       ! aten%u
    real(rk8) , pointer , dimension(:,:,:) :: vten       ! aten%v
    real(rk8) , pointer , dimension(:,:,:,:) :: qxten    ! aten%qx
    real(rk8) , pointer , dimension(:,:,:) :: tketen     ! aten%tke
    real(rk8) , pointer , dimension(:,:,:) :: uuwten     ! uwten%u
    real(rk8) , pointer , dimension(:,:,:) :: vuwten     ! uwten%v
    real(rk8) , pointer , dimension(:,:,:) :: tuwten     ! uwten%t
    real(rk8) , pointer , dimension(:,:,:) :: tkeuwten   ! uwten%tke
    real(rk8) , pointer , dimension(:,:,:,:) :: qxuwten  ! uwten%qx
    real(rk8) , pointer , dimension(:,:,:) :: difft      ! adf%difft
    real(rk8) , pointer , dimension(:,:,:,:) :: diffqx   ! adf%diffqx
    real(rk8) , pointer , dimension(:,:,:,:) :: diagqx   ! holtten%qx
    real(rk8) , pointer , dimension(:,:,:,:) :: chiten   ! chiten
    real(rk8) , pointer , dimension(:,:,:) :: remdrd     ! remdrd
    real(rk8) , pointer , dimension(:,:) :: zpbl
    integer(ik4) , pointer , dimension(:,:) :: kpbl
  end type pbl_2_mod

end module mod_regcm_types

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
