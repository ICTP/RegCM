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
!
! Storage for all the 3d prognostic variables in two
! timesteps and all the 2d variables and constants
!
  type classmsk
    logical , pointer , dimension(:,:) :: landmsk
    logical , pointer , dimension(:,:) :: lakemsk
    logical , pointer , dimension(:,:) :: ocenmsk
    logical , pointer , dimension(:,:,:) :: landmsk1
    logical , pointer , dimension(:,:,:) :: lakemsk1
    logical , pointer , dimension(:,:,:) :: ocenmsk1
  end type classmsk

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
    real(rk8) , pointer , dimension(:,:) :: snowam
    real(rk8) , pointer , dimension(:,:) :: dhlake
  end type domain

  type atmstate
    real(rk8) , pointer , dimension(:,:,:) :: u
    real(rk8) , pointer , dimension(:,:,:) :: v
    real(rk8) , pointer , dimension(:,:,:) :: t
    real(rk8) , pointer , dimension(:,:,:,:) :: qx
    real(rk8) , pointer , dimension(:,:,:) :: tke
  end type atmstate

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
  
  type surfstate
    real(rk8) , pointer , dimension(:,:) :: psa
    real(rk8) , pointer , dimension(:,:) :: psb
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
    real(rk8) , pointer , dimension(:,:,:) :: ubd3d
    real(rk8) , pointer , dimension(:,:,:) :: vbd3d
    real(rk8) , pointer , dimension(:,:,:) :: rhb3d
    real(rk8) , pointer , dimension(:,:,:) :: qsb3d
    real(rk8) , pointer , dimension(:,:,:,:) :: qxb3d
    real(rk8) , pointer , dimension(:,:,:) :: zq
    real(rk8) , pointer , dimension(:,:,:) :: za
    real(rk8) , pointer , dimension(:,:,:) :: dzq
    real(rk8) , pointer , dimension(:,:) :: rhox2d
    real(rk8) , pointer , dimension(:,:,:) :: tkeb3d
    real(rk8) , pointer , dimension(:,:,:,:) :: chib3d
  end type slice

  type diffx
    real(rk8) , pointer , dimension(:,:,:) :: difft
    real(rk8) , pointer , dimension(:,:,:) :: difuu
    real(rk8) , pointer , dimension(:,:,:) :: difuv
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

  type atm_2_lm
    real(rk8) , pointer , dimension(:,:) :: sfps
    real(rk8) , pointer , dimension(:,:) :: tatm
    real(rk8) , pointer , dimension(:,:) :: thatm
    real(rk8) , pointer , dimension(:,:) :: zatm
    real(rk8) , pointer , dimension(:,:) :: qvatm
    real(rk8) , pointer , dimension(:,:) :: uatm
    real(rk8) , pointer , dimension(:,:) :: vatm
    real(rk8) , pointer , dimension(:,:) :: topo
    real(rk8) , pointer , dimension(:,:) :: xlat
    real(rk8) , pointer , dimension(:,:) :: xlon
    real(rk8) , pointer , dimension(:,:) :: lnd
    real(rk8) , pointer , dimension(:,:) :: snow
    real(rk8) , pointer , dimension(:,:) :: rhox
    real(rk8) , pointer , dimension(:,:) :: hpbl
    real(rk8) , pointer , dimension(:,:) :: hfx
    real(rk8) , pointer , dimension(:,:) :: qfx
    real(rk8) , pointer , dimension(:,:) :: tgb
    real(rk8) , pointer , dimension(:,:) :: tgbb
    real(rk8) , pointer , dimension(:,:) :: cosz
    real(rk8) , pointer , dimension(:,:) :: albvgs
    real(rk8) , pointer , dimension(:,:) :: albvgl
    real(rk8) , pointer , dimension(:,:) :: aldirs
    real(rk8) , pointer , dimension(:,:) :: aldifs
    real(rk8) , pointer , dimension(:,:) :: aldirl
    real(rk8) , pointer , dimension(:,:) :: aldifl
    real(rk8) , pointer , dimension(:,:,:) :: topo_sub
    real(rk8) , pointer , dimension(:,:,:) :: lnd_sub
  end type atm_2_lm

  type rcm_2_ocn
    real(rk8) , pointer , dimension(:,:) :: sfps
    real(rk8) , pointer , dimension(:,:) :: tatm
    real(rk8) , pointer , dimension(:,:) :: thatm
    real(rk8) , pointer , dimension(:,:) :: zatm
    real(rk8) , pointer , dimension(:,:) :: qvatm
    real(rk8) , pointer , dimension(:,:) :: uatm
    real(rk8) , pointer , dimension(:,:) :: vatm
    real(rk8) , pointer , dimension(:,:) :: ncprp
    real(rk8) , pointer , dimension(:,:) :: cprp
    real(rk8) , pointer , dimension(:,:) :: rhox
    real(rk8) , pointer , dimension(:,:) :: tgrd
    real(rk8) , pointer , dimension(:,:) :: sst
    logical , pointer , dimension(:,:) :: icemsk
    real(rk8) , pointer , dimension(:,:) :: fsw
    real(rk8) , pointer , dimension(:,:) :: flw
    real(rk8) , pointer , dimension(:,:) :: flwd
    real(rk8) , pointer , dimension(:,:) :: hpbl
    real(rk8) , pointer , dimension(:,:) :: cosz
  end type rcm_2_ocn

  type ocn_2_rcm
    real(rk8) , pointer , dimension(:,:) :: tgrd
    real(rk8) , pointer , dimension(:,:) :: drag
    real(rk8) , pointer , dimension(:,:) :: evpr
    real(rk8) , pointer , dimension(:,:) :: sent
    real(rk8) , pointer , dimension(:,:) :: emiss
    real(rk8) , pointer , dimension(:,:) :: u10m
    real(rk8) , pointer , dimension(:,:) :: v10m
    real(rk8) , pointer , dimension(:,:) :: t2m
    real(rk8) , pointer , dimension(:,:) :: q2m
    real(rk8) , pointer , dimension(:,:) :: tgbrd
    real(rk8) , pointer , dimension(:,:) :: taux
    real(rk8) , pointer , dimension(:,:) :: tauy
    real(rk8) , pointer , dimension(:,:) :: sfice
    real(rk8) , pointer , dimension(:,:) :: sm
    real(rk8) , pointer , dimension(:,:) :: tm
    logical , pointer , dimension(:,:) :: icemsk
    real(rk8) , pointer , dimension(:,:) :: aldirs
    real(rk8) , pointer , dimension(:,:) :: aldifs
    real(rk8) , pointer , dimension(:,:) :: aldirl
    real(rk8) , pointer , dimension(:,:) :: aldifl
  end type ocn_2_rcm

  type ocn_2_sav
    real(rk8) , pointer , dimension(:,:) :: snag
    real(rk8) , pointer , dimension(:,:) :: sncv
    real(rk8) , pointer , dimension(:,:) :: scvk
    real(rk8) , pointer , dimension(:,:) :: deltas
    real(rk8) , pointer , dimension(:,:) :: tdeltas
    real(rk8) , pointer , dimension(:,:) :: dtskin
  end type ocn_2_sav

  type static_2_rad
    real(rk8) , pointer , dimension(:,:) :: xlat
  end type static_2_rad

  type rcm_2_rad
    real(rk8) , pointer , dimension(:,:,:) :: tatms
    real(rk8) , pointer , dimension(:,:,:) :: rhatms
    real(rk8) , pointer , dimension(:,:) :: tground
    real(rk8) , pointer , dimension(:,:) :: psa
    real(rk8) , pointer , dimension(:,:) :: psb
    real(rk8) , pointer , dimension(:,:) :: cosz
    real(rk8) , pointer , dimension(:,:) :: aldirs
    real(rk8) , pointer , dimension(:,:) :: aldifs
    real(rk8) , pointer , dimension(:,:) :: aldirl
    real(rk8) , pointer , dimension(:,:) :: aldifl
    real(rk8) , pointer , dimension(:,:) :: swalb
    real(rk8) , pointer , dimension(:,:) :: lwalb
    real(rk8) , pointer , dimension(:,:) :: emiss
    real(rk8) , pointer , dimension(:,:,:,:) :: qxatms
    real(rk8) , pointer , dimension(:,:,:,:) :: chia
    integer(ik4) , pointer , intent(in) , dimension(:,:) :: ldmsk
  end type rcm_2_rad

  type rad_2_rcm
    real(rk8) , pointer , dimension(:,:) :: abveg
    real(rk8) , pointer , dimension(:,:) :: solar
    real(rk8) , pointer , dimension(:,:) :: totsol
    real(rk8) , pointer , dimension(:,:) :: soldir
    real(rk8) , pointer , dimension(:,:) :: soldif
    real(rk8) , pointer , dimension(:,:) :: netswflx
    real(rk8) , pointer , dimension(:,:) :: netlwflx
    real(rk8) , pointer , dimension(:,:) :: uplwflx
    real(rk8) , pointer , dimension(:,:) :: downlwflx
  end type rad_2_rcm

end module mod_regcm_types
