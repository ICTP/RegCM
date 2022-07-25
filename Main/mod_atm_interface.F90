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

module mod_atm_interface

  use mod_dynparam
  use mod_stdio
  use mod_constants
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_service
  use mod_memutil
  use mod_regcm_types

  implicit none

  private

  logical , public , parameter :: cross = .false.
  logical , public , parameter :: dot = .true.

  type(domain) , public :: mddom
  type(domain_subgrid) , public :: mdsub
  type(atmstate_a) , public :: atm1
  type(atmstate_b) , public :: atm2
  type(atmstate_c) , public :: atmc
  type(atmstate_tendency) , public :: aten
  type(atmstate_decoupled) , public :: atmx
  type(tendiag) , public :: tdiag
  type(qendiag) , public :: qdiag
  type(surfstate) , public :: sfs
  type(slice) , public :: atms
  type(v3dbound) , public :: xtb , xqb , xub , xvb , xppb , xwwb , xpaib
  type(v3dbound) , public :: xlb , xib
  type(v2dbound) , public :: xpsb , xtsb
  type(bound_area) , public :: ba_cr , ba_dt , ba_ut , ba_vt
  type(reference_atmosphere) , public :: atm0
  type(mass_divergence) , public :: mdv
  type(nhboundhelp) , public :: nhbh0 , nhbh1

  ! For idynamic 3

  type(atmosphere) , public :: mo_atm

  public :: allocate_mod_atm_interface
  public :: allocate_v3dbound , allocate_v2dbound
  public :: setup_boundaries , setup_model_indexes
  public :: export_data_from_atm

  real(rkx) , public , pointer , dimension(:,:,:) :: dstor
  real(rkx) , public , pointer , dimension(:,:,:) :: hstor

  real(rkx) , public , pointer , dimension(:,:,:) :: qdot , omega

  ! Sun
  ! Cosine of zenithal solar angle
  real(rkx) , public , pointer , dimension(:,:) :: coszrs

  ! Cumulus
  integer(ik4) , pointer , public , dimension(:,:) :: icumbot
  integer(ik4) , pointer , public , dimension(:,:) :: icumtop
  integer(ik4) , pointer , public , dimension(:,:) :: ktrop
  real(rkx) , pointer , public , dimension(:,:,:) :: convpr
  real(rkx) , pointer , public , dimension(:,:) :: pptc
  real(rkx) , pointer , public , dimension(:,:) :: prca

  ! Radiation
  real(rkx) , pointer , public , dimension(:,:) :: ptrop
  ! vegetation absorbed radiation (full solar spectrum)
  real(rkx) , pointer , public , dimension(:,:) :: sabveg
  ! Incident solar flux
  real(rkx) , pointer , public , dimension(:,:) :: dsol
  real(rkx) , pointer , public , dimension(:,:) :: solis
  real(rkx) , pointer , public , dimension(:,:) :: solvs
  real(rkx) , pointer , public , dimension(:,:) :: solvsd
  real(rkx) , pointer , public , dimension(:,:) :: solvl
  real(rkx) , pointer , public , dimension(:,:) :: solvld
  real(rkx) , pointer , public , dimension(:,:) :: totcf
  real(rkx) , pointer , public , dimension(:,:) :: flw
  real(rkx) , pointer , public , dimension(:,:) :: fsw
  real(rkx) , pointer , public , dimension(:,:) :: flwd
  real(rkx) , pointer , public , dimension(:,:,:) :: cldfra
  real(rkx) , pointer , public , dimension(:,:,:) :: cldlwc
  real(rkx) , pointer , public , dimension(:,:,:) :: heatrt

  ! Dynamic 2
  real(rkx) , pointer , public , dimension(:,:) :: dpsdxm , dpsdym
  real(rkx) , public , dimension(-6:6,-6:6) :: tmask

  ! Surface
  ! Total Long wave albedo (0.7-5.0 micro-meter)
  real(rkx) , pointer , public , dimension(:,:) :: albvl
  ! Total Short wave albedo (0.2-0.7 micro-meter)
  real(rkx) , pointer , public , dimension(:,:) :: albvs
  ! 0.2-0.7 micro-meter srfc alb to direct radiation
  real(rkx) , pointer , public , dimension(:,:) :: aldirs
  ! 0.2-0.7 micro-meter srfc alb to diffuse radiation
  real(rkx) , pointer , public , dimension(:,:) :: aldifs
  ! 0.7-5.0 micro-meter srfc alb to direct radiation
  real(rkx) , pointer , public , dimension(:,:) :: aldirl
  ! 0.7-5.0 micro-meter srfc alb to diffuse radiation
  real(rkx) , pointer , public , dimension(:,:) :: aldifl
  ! Emissivity at surface
  real(rkx) , pointer , public , dimension(:,:) :: emiss
  ! Total solar incoming radiation
  real(rkx) , pointer , public , dimension(:,:) :: sinc

  ! Precip
  real(rkx) , pointer , public , dimension(:,:) :: pptnc
  real(rkx) , pointer , public , dimension(:,:) :: prnca
  real(rkx) , pointer , public , dimension(:,:) :: crrate
  real(rkx) , pointer , public , dimension(:,:) :: ncrrate
  real(rkx) , pointer , public , dimension(:,:,:) :: fcc
  real(rkx) , pointer , public , dimension(:,:,:) :: remrat
  real(rkx) , pointer , public , dimension(:,:,:) :: rembc
  real(rkx) , pointer , public , dimension(:,:,:) :: ccn
  real(rkx) , pointer , public , dimension(:,:,:) :: rain_ls

  ! PBL
  integer(ik4) , public , pointer , dimension(:,:) :: kpbl
  real(rkx) , public , pointer , dimension(:,:) :: zpbl

  ! Cumulus
  real(rkx) , public , pointer , dimension(:,:,:) :: q_detr
  real(rkx) , public , pointer , dimension(:,:,:) :: rain_cc

  ! Surface for chemistry
  real(rkx) , pointer , public , dimension(:,:) :: sdelq
  real(rkx) , pointer , public , dimension(:,:) :: sdelt
  real(rkx) , public , pointer , dimension(:,:) :: ssw2da
  real(rkx) , public , pointer , dimension(:,:) :: sfracv2d
  real(rkx) , public , pointer , dimension(:,:) :: sfracb2d
  real(rkx) , public , pointer , dimension(:,:) :: sfracs2d
  real(rkx) , public , pointer , dimension(:,:) :: svegfrac2d
  real(rkx) , public , pointer , dimension(:,:) :: sxlai2d

#ifdef CLM45
  ! real(rkx) , public , pointer , dimension(:,:) :: ustar
  real(rkx) , public , pointer , dimension(:,:,:) :: voc_em_clm
  real(rkx) , public , pointer , dimension(:,:,:) :: dustflx_clm
  real(rkx) , public , pointer , dimension(:,:,:) :: sw_vol
  real(rkx) , public , pointer , dimension(:,:,:) :: tsoi
#endif

  !chemistry for surface
  real(rkx) , public , pointer , dimension(:,:,:) :: wetdepflx
  real(rkx) , public , pointer , dimension(:,:,:) :: drydepflx

  ! Coupling
  real(rkx) , public , pointer , dimension(:,:,:) :: dailyrnf
  integer(ik4) , public , pointer , dimension(:,:) :: cplmsk

  integer(ik4) :: ix1 , ix2 , jx1 , jx2
  integer(ik4) :: id1 , id2 , jd1 , jd2

#ifdef DEBUG
  !type(grid_nc_var4d) , public , save :: nc_4d
  !type(grid_nc_var3d) , public , save :: nc_3d
  !type(grid_nc_var2d) , public , save :: nc_2d
  !type(grid_nc_var4d) , public , save :: qqxp
#endif

  contains

    subroutine setup_model_indexes
      implicit none
      ma%jbl1 = 1
      ma%jbl2 = 2
      ma%jbl3 = 3
      ma%jbl4 = 4
      ma%jbl6 = 6
      ma%jbr1 = 1
      ma%jbr2 = 2
      ma%jbr3 = 3
      ma%jbr4 = 4
      ma%jbr6 = 6
      ma%ibt1 = 1
      ma%ibt2 = 2
      ma%ibt3 = 3
      ma%ibt4 = 4
      ma%ibt6 = 6
      ma%ibb1 = 1
      ma%ibb2 = 2
      ma%ibb3 = 3
      ma%ibb4 = 4
      ma%ibb6 = 6
      if ( ma%has_bdyleft ) then
        ma%jbl1 = 0
        ma%jbl2 = 0
        ma%jbl3 = 0
        ma%jbl4 = 0
        ma%jbl6 = 0
      end if
      if ( ma%has_bdyright ) then
        ma%jbr1 = 0
        ma%jbr2 = 0
        ma%jbr3 = 0
        ma%jbr4 = 0
        ma%jbr6 = 0
      end if
      if ( ma%has_bdytop ) then
        ma%ibt1 = 0
        ma%ibt2 = 0
        ma%ibt3 = 0
        ma%ibt4 = 0
        ma%ibt6 = 0
      end if
      if ( ma%has_bdybottom ) then
        ma%ibb1 = 0
        ma%ibb2 = 0
        ma%ibb3 = 0
        ma%ibb4 = 0
        ma%ibb6 = 0
      end if
      jde1  = global_dot_jstart
      jdi1  = global_dot_jstart
      jdii1 = global_dot_jstart
      jde2  = global_dot_jend
      jdi2  = global_dot_jend
      jdii2 = global_dot_jend
      ide1  = global_dot_istart
      idi1  = global_dot_istart
      idii1 = global_dot_istart
      ide2  = global_dot_iend
      idi2  = global_dot_iend
      idii2 = global_dot_iend
      if ( ma%has_bdyleft ) then
        jdi1 = jde1 + 1
        jdii1 = jde1 + 2
      end if
      if ( ma%has_bdyright ) then
        jdi2 = jde2 - 1
        jdii2 = jde2 - 2
      end if
      if ( ma%has_bdybottom ) then
        idi1 = ide1 + 1
        idii1 = ide1 + 2
      end if
      if ( ma%has_bdytop ) then
        idi2 = ide2 - 1
        idii2 = ide2 - 2
      end if
      jce1  = global_cross_jstart
      jci1  = global_cross_jstart
      jcii1 = global_cross_jstart
      jce2  = global_cross_jend
      jci2  = global_cross_jend
      jcii2 = global_cross_jend
      ice1  = global_cross_istart
      ici1  = global_cross_istart
      icii1 = global_cross_istart
      ice2  = global_cross_iend
      ici2  = global_cross_iend
      icii2 = global_cross_iend
      if ( ma%has_bdyleft ) then
        jci1 = jce1 + 1
        jcii1 = jce1 + 2
      end if
      if ( ma%has_bdyright ) then
        jci2 = jce2 - 1
        jcii2 = jce2 - 2
      end if
      if ( ma%has_bdybottom ) then
        ici1 = ice1 + 1
        icii1 = ice1 + 2
      end if
      if ( ma%has_bdytop ) then
        ici2 = ice2 - 1
        icii2 = ice2 - 2
      end if
      idi1ga = idi1 - ma%ibb1
      idi2ga = idi2 + ma%ibt1
      jdi1ga = jdi1 - ma%jbl1
      jdi2ga = jdi2 + ma%jbr1
      ide1ga = ide1 - ma%ibb1
      ide2ga = ide2 + ma%ibt1
      jde1ga = jde1 - ma%jbl1
      jde2ga = jde2 + ma%jbr1
      ici1ga = ici1 - ma%ibb1
      ici2ga = ici2 + ma%ibt1
      jci1ga = jci1 - ma%jbl1
      jci2ga = jci2 + ma%jbr1
      ice1ga = ice1 - ma%ibb1
      ice2ga = ice2 + ma%ibt1
      jce1ga = jce1 - ma%jbl1
      jce2ga = jce2 + ma%jbr1

      idi1gb = idi1 - ma%ibb2
      idi2gb = idi2 + ma%ibt2
      jdi1gb = jdi1 - ma%jbl2
      jdi2gb = jdi2 + ma%jbr2
      ide1gb = ide1 - ma%ibb2
      ide2gb = ide2 + ma%ibt2
      jde1gb = jde1 - ma%jbl2
      jde2gb = jde2 + ma%jbr2
      ici1gb = ici1 - ma%ibb2
      ici2gb = ici2 + ma%ibt2
      jci1gb = jci1 - ma%jbl2
      jci2gb = jci2 + ma%jbr2
      ice1gb = ice1 - ma%ibb2
      ice2gb = ice2 + ma%ibt2
      jce1gb = jce1 - ma%jbl2
      jce2gb = jce2 + ma%jbr2

      idi1gc = idi1 - ma%ibb3
      idi2gc = idi2 + ma%ibt3
      jdi1gc = jdi1 - ma%jbl3
      jdi2gc = jdi2 + ma%jbr3
      ide1gc = ide1 - ma%ibb3
      ide2gc = ide2 + ma%ibt3
      jde1gc = jde1 - ma%jbl3
      jde2gc = jde2 + ma%jbr3
      ici1gc = ici1 - ma%ibb3
      ici2gc = ici2 + ma%ibt3
      jci1gc = jci1 - ma%jbl3
      jci2gc = jci2 + ma%jbr3
      ice1gc = ice1 - ma%ibb3
      ice2gc = ice2 + ma%ibt3
      jce1gc = jce1 - ma%jbl3
      jce2gc = jce2 + ma%jbr3

      idi1sl = idi1gb - ma%ibb2
      idi2sl = idi2gb + ma%ibt2
      jdi1sl = jdi1gb - ma%jbl2
      jdi2sl = jdi2gb + ma%jbr2
      ide1sl = ide1gb - ma%ibb2
      ide2sl = ide2gb + ma%ibt2
      jde1sl = jde1gb - ma%jbl2
      jde2sl = jde2gb + ma%jbr2
      ici1sl = ici1gb - ma%ibb2
      ici2sl = ici2gb + ma%ibt2
      jci1sl = jci1gb - ma%jbl2
      jci2sl = jci2gb + ma%jbr2
      ice1sl = ice1gb - ma%ibb2
      ice2sl = ice2gb + ma%ibt2
      jce1sl = jce1gb - ma%jbl2
      jce2sl = jce2gb + ma%jbr2

      jde1sg = (jde1-1)*nsg+1
      jde2sg = jde2*nsg
      ide1sg = (ide1-1)*nsg+1
      ide2sg = ide2*nsg

#ifdef DEBUG
      write(ndebug,*) 'TOPLEFT     = ', ma%topleft
      write(ndebug,*) 'TOP         = ', ma%top
      write(ndebug,*) 'TOPRIGHT    = ', ma%topright
      write(ndebug,*) 'RIGHT       = ', ma%right
      write(ndebug,*) 'BOTTOMRIGHT = ', ma%bottomright
      write(ndebug,*) 'BOTTOM      = ', ma%bottom
      write(ndebug,*) 'BOTTOMLEFT  = ', ma%bottomleft
      write(ndebug,*) 'LEFT        = ', ma%left
      write(ndebug,*) 'DOTPEXTJI1 : ', jde1 , jde2 , ide1 , ide2
      write(ndebug,*) 'DOTPINTJI1 : ', jdi1 , jdi2 , idi1 , idi2
      write(ndebug,*) 'DOTPINTJI2 : ', jdii1 , jdii2 , idii1 , idii2
      write(ndebug,*) 'CRXPEXTJI1 : ', jce1 , jce2 , ice1 , ice2
      write(ndebug,*) 'CRXPINTJI1 : ', jci1 , jci2 , ici1 , ici2
      write(ndebug,*) 'CRXPINTJI2 : ', jcii1 , jcii2 , icii1 , icii2
      write(ndebug,*) 'TOPBDY   : ', ma%has_bdytop
      write(ndebug,*) 'BTMBDY   : ', ma%has_bdybottom
      write(ndebug,*) 'RGTBDY   : ', ma%has_bdyright
      write(ndebug,*) 'LFTBDY   : ', ma%has_bdyleft
      flush(ndebug)
#endif
    end subroutine setup_model_indexes

    subroutine setup_boundaries(ldotx,ldoty,ba)
      implicit none
      logical , intent(in) :: ldotx , ldoty
      type(bound_area) , intent(out) :: ba
      integer(ik4) :: icx , icy
      integer(ik4) :: igbb1 , igbb2 , igbt1 , igbt2
      integer(ik4) :: jgbl1 , jgbl2 , jgbr1 , jgbr2
      integer(ik4) :: i , j , i1 , i2 , j1 , j2

      call getmem2d(ba%ibnd,jde1,jde2,ide1,ide2,'setup_boundaries:ibnd')
      call getmem2d(ba%bsouth,jde1,jde2,ide1,ide2,'setup_boundaries:bsouth')
      call getmem2d(ba%bnorth,jde1,jde2,ide1,ide2,'setup_boundaries:bnorth')
      call getmem2d(ba%beast,jde1,jde2,ide1,ide2,'setup_boundaries:beast')
      call getmem2d(ba%bwest,jde1,jde2,ide1,ide2,'setup_boundaries:bwest')

      if ( ldotx ) then
        icx = 0
      else
        icx = 1
      end if
      if ( ldoty ) then
        icy = 0
      else
        icy = 1
      end if
      if ( ldotx .or. ldoty ) then
        ba%nsp = nspgd
      else
        ba%nsp = nspgx
      end if
      ba%ibnd(:,:) = -1
      igbb1 = 2
      igbb2 = ba%nsp-1
      jgbl1 = 2
      jgbl2 = ba%nsp-1
      igbt1 = iy-icy-ba%nsp+2
      igbt2 = iy-icy-1
      jgbr1 = jx-icx-ba%nsp+2
      jgbr2 = jx-icx-1
      i1 = ide1
      i2 = ide2
      j1 = jde1
      j2 = jde2
      ba%ns = 0
      ba%nn = 0
      ba%nw = 0
      ba%ne = 0
      ba%bsouth(:,:) = .false.
      ba%bnorth(:,:) = .false.
      ba%bwest(:,:) = .false.
      ba%beast(:,:) = .false.
      if (.not. ma%crmflag ) then
        if ( ma%bandflag ) then
          ! Check for South boundary
          do i = i1 , i2
            if ( i >= igbb1 .and. i <= igbb2 ) then
              do j = j1 , j2
                if ( j < jgbl1 .and. j > jgbr2 ) cycle
                ba%ibnd(j,i) = i-igbb1+2
                ba%bsouth(j,i) = .true.
                ba%ns = ba%ns+1
              end do
            end if
          end do
          ! North Boundary
          do i = i1 , i2
            if ( i >= igbt1 .and. i <= igbt2 ) then
              do j = j1 , j2
                if ( j < jgbl1 .and. j > jgbr2 ) cycle
                ba%ibnd(j,i) = igbt2-i+2
                ba%bnorth(j,i) = .true.
                ba%nn = ba%nn+1
              end do
            end if
          end do
        else
          ! Check for South boundary
          do i = i1 , i2
            if ( i >= igbb1 .and. i <= igbb2 ) then
              do j = j1 , j2
                if ( j >= jgbl1 .and. j <= jgbr2 ) then
                  if ( j <= jgbl2 .and. i >= j ) cycle
                  if ( j >= jgbr1 .and. i >= (jgbr2-j+2) ) cycle
                  ba%ibnd(j,i) = i-igbb1+2
                  ba%bsouth(j,i) = .true.
                  ba%ns = ba%ns+1
                end if
              end do
            end if
          end do
          ! North Boundary
          do i = i1 , i2
            if ( i >= igbt1 .and. i <= igbt2 ) then
              do j = j1 , j2
                if ( j >= jgbl1 .and. j <= jgbr2 ) then
                  if ( j <= jgbl2 .and. (igbt2-i+2) >= j ) cycle
                  if ( j >= jgbr1 .and. (igbt2-i) >= (jgbr2-j) ) cycle
                  ba%ibnd(j,i) = igbt2-i+2
                  ba%bnorth(j,i) = .true.
                  ba%nn = ba%nn+1
                end if
              end do
            end if
          end do
          ! West boundary
          do i = i1 , i2
            if ( i < igbb1 .or. i > igbt2 ) cycle
            do j = j1 , j2
              if ( j >= jgbl1 .and. j <= jgbl2 ) then
                if ( i < igbb2 .and. j > i ) cycle
                if ( i > igbt1 .and. j > (igbt2-i+2) ) cycle
                ba%ibnd(j,i) = j-jgbl1+2
                ba%bwest(j,i) = .true.
                ba%nw = ba%nw+1
              end if
            end do
          end do
          ! East boundary
          do i = i1 , i2
            if ( i < igbb1 .or. i > igbt2 ) cycle
            do j = j1 , j2
              if ( j >= jgbr1 .and. j <= jgbr2 ) then
                if ( i < igbb2 .and. (jgbr2-j+2) > i ) cycle
                if ( i > igbt1 .and. (jgbr2-j) > (igbt2-i) ) cycle
                ba%ibnd(j,i) = jgbr2-j+2
                ba%beast(j,i) = .true.
                ba%ne = ba%ne+1
              end if
            end do
          end do
        end if
      end if
      ba%havebound = (ba%ns /= 0 .or. ba%nn /= 0 .or. &
                      ba%nw /= 0 .or. ba%ne /= 0)
#ifdef DEBUG
      write(ndebug,*) 'DOTX : ', ldotx
      write(ndebug,*) 'DOTY : ', ldoty
      write(ndebug,*) 'BDYS : ', ba%ns
      write(ndebug,*) 'BDYN : ', ba%nn
      write(ndebug,*) 'BDYW : ', ba%nw
      write(ndebug,*) 'BDYE : ', ba%ne

      do i = ide2 , ide1 , -1
        do j = jde1, jde2
          if ( ba%bsouth(j,i) ) then
            write(ndebug,'(1a,i0.4)',advance='no') 'S' , ba%ibnd(j,i)
          else if ( ba%bnorth(j,i) ) then
            write(ndebug,'(1a,i0.4)',advance='no') 'N' , ba%ibnd(j,i)
          else if ( ba%bwest(j,i) ) then
            write(ndebug,'(1a,i0.4)',advance='no') 'W' , ba%ibnd(j,i)
          else if ( ba%beast(j,i) ) then
            write(ndebug,'(1a,i0.4)',advance='no') 'E' , ba%ibnd(j,i)
          else
            write(ndebug,'(1a,i0.4)',advance='no') 'X', 0
          end if
        end do
        write(ndebug,*) ' '
      end do
#endif
    end subroutine setup_boundaries

    subroutine allocate_v3dbound(xb,ke,ldot)
      implicit none
      type(v3dbound) , intent(out) :: xb
      integer(ik4) , intent(in) :: ke
      logical , intent(in) :: ldot
      if ( ldot ) then
        call getmem3d(xb%b0,jde1ga,jde2ga,ide1ga,ide2ga,1,ke,'v3dbound:b0')
        call getmem3d(xb%b1,jde1ga,jde2ga,ide1ga,ide2ga,1,ke,'v3dbound:b1')
        call getmem3d(xb%bt,jde1ga,jde2ga,ide1ga,ide2ga,1,ke,'v3dbound:bt')
      else
        call getmem3d(xb%b0,jce1ga,jce2ga,ice1ga,ice2ga,1,ke,'v3dbound:b0')
        call getmem3d(xb%b1,jce1ga,jce2ga,ice1ga,ice2ga,1,ke,'v3dbound:b1')
        call getmem3d(xb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,ke,'v3dbound:bt')
      end if
    end subroutine allocate_v3dbound

    subroutine allocate_v2dbound(xb,ldot)
      implicit none
      type(v2dbound) , intent(out) :: xb
      logical , intent(in) :: ldot
      if ( ldot ) then
        call getmem2d(xb%b0,jde1ga,jde2ga,ide1ga,ide2ga,'v2dbound:b0')
        call getmem2d(xb%b1,jde1ga,jde2ga,ide1ga,ide2ga,'v2dbound:b1')
        call getmem2d(xb%bt,jde1ga,jde2ga,ide1ga,ide2ga,'v2dbound:bt')
      else
        call getmem2d(xb%b0,jce1ga,jce2ga,ice1ga,ice2ga,'v2dbound:b0')
        call getmem2d(xb%b1,jce1ga,jce2ga,ice1ga,ice2ga,'v2dbound:b1')
        call getmem2d(xb%bt,jce1ga,jce2ga,ice1ga,ice2ga,'v2dbound:bt')
      end if
    end subroutine allocate_v2dbound

    subroutine allocate_atmosphere(atm)
      implicit none
      type(atmosphere) , intent(inout) :: atm
      call getmem3d(atm%u,jde1gb,jde2gb,ice1ga,ice2ga,1,kz,'atmstate:u')
      call getmem3d(atm%v,jce1ga,jce2ga,ide1gb,ide2gb,1,kz,'atmstate:v')
      call getmem3d(atm%ux,jce1gb,jce2gb,ice1ga,ice2ga,1,kz,'atmstate:ux')
      call getmem3d(atm%vx,jce1ga,jce2ga,ice1gb,ice2gb,1,kz,'atmstate:vx')
      call getmem3d(atm%w,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:w')
      call getmem3d(atm%pai,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:pai')
      call getmem3d(atm%p,jce1,jce2,ice1,ice2,1,kz,'atmstate:p')
      call getmem3d(atm%rho,jce1,jce2,ice1,ice2,1,kz,'atmstate:rho')
      call getmem3d(atm%pf,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:pf')
      call getmem3d(atm%t,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:t')
      call getmem3d(atm%tvirt,jce1,jce2,ice1,ice2,1,kz,'atmstate:tvirt')
      call getmem3d(atm%tetav,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:tetav')
      call getmem3d(atm%zeta,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'atmstate:zeta')
      call getmem3d(atm%zetaf,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:zetaf')
      call getmem3d(atm%dz,jce1,jce2,ice1,ice2,1,kz,'atmstate:dz')
      call getmem3d(atm%qs,jce1,jce2,ice1,ice2,1,kz,'atmstate:qs')
      call getmem4d(atm%qx,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,1,nqx,'atmstate:qx')
      call getmem3d(atm%tten,jci1,jci2,ici1,ici2,1,kz,'atmstate:tten')
      call getmem3d(atm%uten,jdi1,jdi2,ici1,ici2,1,kz,'atmstate:uten')
      call getmem3d(atm%vten,jci1,jci2,idi1,idi2,1,kz,'atmstate:vten')
      call getmem4d(atm%qxten,jci1,jci2,ici1,ici2,1,kz,1,nqx,'atmstate:qxten')
      if ( ibltyp == 2 ) then
        call getmem3d(atm%tke,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:tke')
        call getmem3d(atm%tketen,jci1,jci2,ici1,ici2,1,kzp1,'atmstate:tketen')
      end if
      if ( ichem == 1 ) then
        call getmem4d(atm%trac,jce1ga,jce2ga, &
                               ice1ga,ice2ga,1,kz,1,ntr,'atmstate:trac')
        call getmem4d(atm%chiten,jci1,jci2,   &
                                 ici1,ici2,1,kz,1,ntr,'atmstate:chiten')
      end if
      call getmem3d(atm%fmz,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:fmz')
      call getmem3d(atm%fmzf,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:fmzf')
    end subroutine allocate_atmosphere

    subroutine allocate_atmstate_a(atm)
      implicit none
      type(atmstate_a) , intent(out) :: atm
      call getmem3d(atm%u,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:u')
      call getmem3d(atm%v,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:v')
      call getmem3d(atm%t,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:t')
      call getmem4d(atm%qx,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,1,nqx,'atmstate:qx')
      if ( ibltyp == 2 ) then
        call getmem3d(atm%tke,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem3d(atm%pr,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:rho')
        call getmem3d(atm%pp,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:pp')
        call getmem3d(atm%w,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'atmstate:w')
      else
        call getmem3d(atm%pr,jce1,jce2,ice1,ice2,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1,jce2,ice1,ice2,1,kz,'atmstate:rho')
      end if
      if ( ichem == 1 ) then
        call getmem4d(atm%chi,jce1ga,jce2ga,ice1ga,ice2ga, &
                              1,kz,1,ntr,'atmstate:chi')
      end if
    end subroutine allocate_atmstate_a

    subroutine allocate_atmstate_b(atm)
      implicit none
      type(atmstate_b) , intent(out) :: atm

      call getmem3d(atm%u,jd1,jd2,id1,id2,1,kz,'atmstate:u')
      call getmem3d(atm%v,jd1,jd2,id1,id2,1,kz,'atmstate:v')
      call getmem3d(atm%t,jx1,jx2,ix1,ix2,1,kz,'atmstate:t')
      if ( isladvec == 1 ) then
        call getmem4d(atm%qx,jce1sl,jce2sl, &
                             ice1sl,ice2sl,1,kz,1,nqx,'atmstate:qx')
      else
        call getmem4d(atm%qx,jx1,jx2,ix1,ix2,1,kz,1,nqx,'atmstate:qx')
      end if
      if ( ibltyp == 2 ) then
        call getmem3d(atm%tke,jx1,jx2,ix1,ix2,1,kzp1,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem3d(atm%pp,jx1,jx2,ix1,ix2,1,kz,'atmstate:pp')
        call getmem3d(atm%w,jx1,jx2,ix1,ix2,1,kzp1,'atmstate:w')
      end if
      if ( ichem == 1 ) then
        if ( isladvec == 1 ) then
          call getmem4d(atm%chi,jce1sl,jce2sl,ice1sl,ice2sl, &
                                1,kz,1,ntr,'atmstate:chi')
        else
          call getmem4d(atm%chi,jx1,jx2,ix1,ix2,1,kz,1,ntr,'atmstate:chi')
        end if
      end if
      call getmem3d(atm%pr,jce1,jce2,ice1,ice2,1,kz,'atmstate:pr')
    end subroutine allocate_atmstate_b

    subroutine allocate_atmstate_c(atm)
      implicit none
      type(atmstate_c) , intent(out) :: atm
      if ( ibltyp == 2 ) then
        call getmem3d(atm%tke,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem3d(atm%u,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:u')
        call getmem3d(atm%v,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:v')
        call getmem3d(atm%pp,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:pp')
        call getmem3d(atm%t,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:t')
        call getmem3d(atm%w,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:w')
      else
        call getmem3d(atm%u,jdi1,jdi2,idi1,idi2,1,kz,'atmstate:u')
        call getmem3d(atm%v,jdi1,jdi2,idi1,idi2,1,kz,'atmstate:v')
        call getmem3d(atm%t,jci1,jci2,ici1,ici2,1,kz,'atmstate:t')
      end if
      call getmem4d(atm%qx,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,1,nqx,'atmstate:qx')
      if ( ichem == 1 ) then
        call getmem4d(atm%chi,jce1ga,jce2ga,ice1ga,ice2ga, &
                              1,kz,1,ntr,'atmstate:chi')
      end if
    end subroutine allocate_atmstate_c

    subroutine allocate_atmstate_decoupled(atm)
      implicit none
      type(atmstate_decoupled) , intent(out) :: atm
      call getmem3d(atm%uc,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:uc')
      call getmem3d(atm%vc,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:vc')
      call getmem3d(atm%umc,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:umc')
      call getmem3d(atm%vmc,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:vmc')
      if ( isladvec == 1 ) then
        call getmem3d(atm%ud,jde1gb,jde2gb,ide1gb,ide2gb,1,kz,'atmstate:ud')
        call getmem3d(atm%vd,jde1gb,jde2gb,ide1gb,ide2gb,1,kz,'atmstate:vd')
        call getmem3d(atm%umd,jde1gb,jde2gb,ide1gb,ide2gb,1,kz,'atmstate:umd')
        call getmem3d(atm%vmd,jde1gb,jde2gb,ide1gb,ide2gb,1,kz,'atmstate:vmd')
      else
        call getmem3d(atm%ud,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:ud')
        call getmem3d(atm%vd,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:vd')
        call getmem3d(atm%umd,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:umd')
        call getmem3d(atm%vmd,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:vmd')
      end if
      call getmem3d(atm%t,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:t')
      call getmem3d(atm%tv,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:tv')
      call getmem4d(atm%qx,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,1,nqx,'atmstate:qx')
      if ( idynamic == 2 ) then
        call getmem3d(atm%pr,jce1,jce2,ice1,ice2,1,kz,'atmstate:pr')
        call getmem3d(atm%pp,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:pp')
        call getmem3d(atm%w,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'atmstate:w')
      end if
      if ( ichem == 1 ) then
        call getmem4d(atm%chi,jce1ga,jce2ga,ice1ga,ice2ga, &
                              1,kz,1,ntr,'atmstate:chi')
      end if
    end subroutine allocate_atmstate_decoupled

    subroutine allocate_atmstate_tendency(atm)
      implicit none
      type(atmstate_tendency) , intent(out) :: atm
      call getmem4d(atm%u,jdi1,jdi2,idi1,idi2, &
                    1,kz,1,number_of_prognostic_components,'atmstate:u')
      call getmem4d(atm%v,jdi1,jdi2,idi1,idi2, &
                    1,kz,1,number_of_prognostic_components,'atmstate:v')
      call getmem4d(atm%t,jci1,jci2,ici1,ici2, &
                    1,kz,1,number_of_prognostic_components,'atmstate:t')
      call getmem5d(atm%qx,jci1,jci2,ici1,ici2, &
                    1,kz,1,nqx,1,number_of_prognostic_components,'atmstate:qx')
      if ( ibltyp == 2 ) then
        call getmem4d(atm%tke,jci1,jci2,ici1,ici2, &
                  1,kzp1,1,number_of_prognostic_components,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem4d(atm%pp,jci1,jci2,ici1,ici2, &
                      1,kz,1,number_of_prognostic_components,'atmstate:pp')
        call getmem4d(atm%w,jci1,jci2,ici1,ici2, &
                      1,kzp1,1,number_of_prognostic_components,'atmstate:w')
      end if
      if ( ichem == 1 ) then
        call getmem5d(atm%chi,jci1,jci2,ici1,ici2, &
                      1,kz,1,ntr,1,number_of_prognostic_components, &
                      'atmstate:chi')
      end if
    end subroutine allocate_atmstate_tendency

    subroutine allocate_reference_atmosphere(atm)
      implicit none
      type(reference_atmosphere) , intent(out) :: atm
      call getmem2d(atm%ps,jce1ga,jce2ga,ice1ga,ice2ga,'reference:ps')
      call getmem2d(atm%psdot,jde1ga,jde2ga,ide1ga,ide2ga,'reference:psdot')
      call getmem3d(atm%t,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'reference:t')
      call getmem3d(atm%pr,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'reference:pr')
      call getmem3d(atm%rho,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'reference:rho')
      call getmem3d(atm%z,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'reference:z')
      call getmem3d(atm%zd,jdi1,jdi2,idi1,idi2,1,kz,'reference:zd')
      call getmem3d(atm%tf,jce1,jce2,ice1,ice2,1,kzp1,'reference:tf')
      call getmem3d(atm%pf,jce1,jce2,ice1,ice2,1,kzp1,'reference:pf')
      call getmem3d(atm%rhof,jce1,jce2,ice1,ice2,1,kzp1,'reference:rhof')
      call getmem3d(atm%zf,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'reference:zf')
      call getmem3d(atm%dzf,jce1,jce2,ice1,ice2,1,kzp1,'reference:dzf')
      call getmem3d(atm%dprddx,jdi1,jdi2,idi1,idi2,1,kz,'reference:dprddx')
      call getmem3d(atm%dprddy,jdi1,jdi2,idi1,idi2,1,kz,'reference:dprddy')
    end subroutine allocate_reference_atmosphere

    subroutine allocate_mass_divergence(div)
      implicit none
      type(mass_divergence) , intent(out) :: div
      call getmem3d(div%cr,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'massdiv:cr')
    end subroutine allocate_mass_divergence

    subroutine allocate_tendiag(dia)
      implicit none
      type(tendiag) , intent(out) :: dia
      call getmem3d(dia%adh,jci1,jci2,ici1,ici2,1,kz,'tendiag:adh')
      call getmem3d(dia%adv,jci1,jci2,ici1,ici2,1,kz,'tendiag:adv')
      call getmem3d(dia%tbl,jci1,jci2,ici1,ici2,1,kz,'tendiag:tbl')
      call getmem3d(dia%con,jci1,jci2,ici1,ici2,1,kz,'tendiag:con')
      call getmem3d(dia%bdy,jci1,jci2,ici1,ici2,1,kz,'tendiag:bdy')
      call getmem3d(dia%adi,jci1,jci2,ici1,ici2,1,kz,'tendiag:adi')
      call getmem3d(dia%dif,jci1,jci2,ici1,ici2,1,kz,'tendiag:dif')
      call getmem3d(dia%rad,jci1,jci2,ici1,ici2,1,kz,'tendiag:rad')
      call getmem3d(dia%lsc,jci1,jci2,ici1,ici2,1,kz,'tendiag:lsc')
    end subroutine allocate_tendiag

    subroutine allocate_qendiag(dia)
      implicit none
      type(qendiag) , intent(out) :: dia
      call getmem3d(dia%adh,jci1,jci2,ici1,ici2,1,kz,'tendiag:adh')
      call getmem3d(dia%adv,jci1,jci2,ici1,ici2,1,kz,'tendiag:adv')
      call getmem3d(dia%tbl,jci1,jci2,ici1,ici2,1,kz,'tendiag:tbl')
      call getmem3d(dia%con,jci1,jci2,ici1,ici2,1,kz,'tendiag:con')
      call getmem3d(dia%bdy,jci1,jci2,ici1,ici2,1,kz,'tendiag:bdy')
      call getmem3d(dia%adi,jci1,jci2,ici1,ici2,1,kz,'tendiag:adi')
      call getmem3d(dia%dif,jci1,jci2,ici1,ici2,1,kz,'tendiag:dif')
      call getmem3d(dia%rad,jci1,jci2,ici1,ici2,1,kz,'tendiag:rad')
      call getmem3d(dia%lsc,jci1,jci2,ici1,ici2,1,kz,'tendiag:lsc')
      if ( ichem == 1 .and. iaerosol == 1 .and. iindirect == 2 ) then
        call getmem3d(dia%qcl,jci1,jci2,ici1,ici2,1,kz,'tendiag:qcl')
        call getmem3d(dia%qcr,jci1,jci2,ici1,ici2,1,kz,'tendiag:qcr')
        call getmem3d(dia%acr,jci1,jci2,ici1,ici2,1,kz,'tendiag:acr')
      end if
    end subroutine allocate_qendiag

    subroutine allocate_domain(dom)
      implicit none
      type(domain) , intent(out) :: dom
      call getmem2d(dom%ht,jde1gb,jde2gb,ide1gb,ide2gb,'storage:ht')
      call getmem2d(dom%lndcat,jde1,jde2,ide1,ide2,'storage:lndcat')
      call getmem2d(dom%lndtex,jde1,jde2,ide1,ide2,'storage:lndtex')
      call getmem2d(dom%xlat,jde1ga,jde2ga,ide1ga,ide2ga,'storage:xlat')
      call getmem2d(dom%xlon,jde1ga,jde2ga,ide1ga,ide2ga,'storage:xlon')
      call getmem2d(dom%dlat,jde1,jde2,ide1,ide2,'storage:dlat')
      call getmem2d(dom%dlon,jde1,jde2,ide1,ide2,'storage:dlon')
      call getmem2d(dom%mask,jde1,jde2,ide1,ide2,'storage:mask')
      call getmem2d(dom%area,jde1,jde2,ide1,ide2,'storage:area')
      if ( idynamic == 3 ) then
        call getmem2d(dom%msfx,jde1,jde2,ide1,ide2,'storage:msfx')
        call getmem2d(dom%msfu,jde1ga,jde2ga,ide1,ide2,'storage:msfu')
        call getmem2d(dom%msfv,jde1,jde2,ide1ga,ide2ga,'storage:msfv')
        call getmem2d(dom%hx,jde1ga,jde2ga,ice1,ice2,'storage:hx')
        call getmem2d(dom%hy,jce1,jce2,ide1ga,ide2ga,'storage:hy')
        call getmem2d(dom%ulat,jde1,jde2,ide1,ide2,'storage:ulat')
        call getmem2d(dom%ulon,jde1,jde2,ide1,ide2,'storage:ulon')
        call getmem2d(dom%vlat,jde1,jde2,ide1,ide2,'storage:vlat')
        call getmem2d(dom%vlon,jde1,jde2,ide1,ide2,'storage:vlon')
        call getmem2d(dom%coriou,jde1,jde2,ice1,ice2,'storage:fu')
        call getmem2d(dom%coriov,jce1,jce2,ide1,ide2,'storage:fv')
      else
        call getmem2d(dom%msfx,jd1,jd2,id1,id2,'storage:msfx')
        call getmem2d(dom%msfd,jd1,jd2,id1,id2,'storage:msfd')
      end if
      call getmem2d(dom%coriol,jde1,jde2,ide1,ide2,'storage:f')
      call getmem2d(dom%snowam,jde1,jde2,ide1,ide2,'storage:snowam')
      call getmem2d(dom%smoist,jde1,jde2,ide1,ide2,'storage:smoist')
      call getmem3d(dom%rmoist,jde1,jde2,ide1,ide2, &
                    1,num_soil_layers,'storage:rmoist')
      call getmem2d(dom%ldmsk,jci1,jci2,ici1,ici2,'storage:ldmsk')
      call getmem2d(dom%iveg,jci1,jci2,ici1,ici2,'storage:iveg')
      call getmem2d(dom%itex,jci1,jci2,ici1,ici2,'storage:itex')
      call getmem2d(dom%xmsf,jdi1,jdi2,idi1,idi2,'storage:xmsf')
      call getmem2d(dom%dmsf,jdi1,jdi2,idi1,idi2,'storage:dmsf')
      if ( lakemod == 1 ) then
        call getmem2d(dom%dhlake,jde1,jde2,ide1,ide2,'storage:dhlake')
      end if
      if ( idynamic == 2 ) then
        call getmem2d(dom%ef,jdi1ga,jdi2ga,idi1ga,idi2ga,'storage:ef')
        call getmem2d(dom%ddx,jdi1ga,jdi2ga,idi1ga,idi2ga,'storage:ddx')
        call getmem2d(dom%ddy,jdi1ga,jdi2ga,idi1ga,idi2ga,'storage:ddy')
        call getmem2d(dom%ex,jci1,jci2,ici1,ici2,'storage:ex')
        call getmem2d(dom%crx,jci1,jci2,ici1,ici2,'storage:crx')
        call getmem2d(dom%cry,jci1,jci2,ici1,ici2,'storage:cry')
        call getmem2d(dom%dmdy,jdi1,jdi2,idi1,idi2,'storage:dmdy')
        call getmem2d(dom%dmdx,jdi1,jdi2,idi1,idi2,'storage:dmdx')
      end if
    end subroutine allocate_domain

    subroutine allocate_domain_subgrid(sub)
      implicit none
      type(domain_subgrid) , intent(out) :: sub
      call getmem3d(sub%ht,1,nnsg,jde1,jde2,ide1,ide2,'storage:ht')
      call getmem3d(sub%lndcat,1,nnsg,jde1,jde2,ide1,ide2,'storage:lndcat')
      call getmem3d(sub%lndtex,1,nnsg,jde1,jde2,ide1,ide2,'storage:lndtex')
      call getmem3d(sub%xlat,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlat')
      call getmem3d(sub%xlon,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlon')
      call getmem3d(sub%mask,1,nnsg,jde1,jde2,ide1,ide2,'storage:mask')
      call getmem3d(sub%area,1,nnsg,jde1,jde2,ide1,ide2,'storage:area')
      call getmem3d(sub%ldmsk,1,nnsg,jci1,jci2,ici1,ici2,'storage:ldmsk')
      call getmem3d(sub%iveg,1,nnsg,jci1,jci2,ici1,ici2,'storage:iveg')
      call getmem3d(sub%itex,1,nnsg,jci1,jci2,ici1,ici2,'storage:itex')
      if ( lakemod == 1 ) then
        call getmem3d(sub%dhlake,1,nnsg,jde1,jde2,ide1,ide2,'storage:dhlake')
      end if
    end subroutine allocate_domain_subgrid

    subroutine allocate_surfstate(sfs)
      implicit none
      type(surfstate) , intent(out) :: sfs
      call getmem2d(sfs%psa,jce1ga,jce2ga,ice1ga,ice2ga,'surf:psa')
      if ( idynamic == 3 ) then
        call assignpnt(sfs%psa,sfs%psb)
        call assignpnt(sfs%psa,sfs%psc)
      else
        call getmem2d(sfs%psdota,jde1ga,jde2ga,ide1ga,ide2ga,'surf:psdota')
        call getmem2d(sfs%psb,jx1,jx2,ix1,ix2,'surf:psb')
        call getmem2d(sfs%psdotb,jd1,jd2,id1,id2,'surf:psdotb')
        call getmem2d(sfs%psc,jce1,jce2,ice1,ice2,'surf:psc')
      end if
      call getmem2d(sfs%tg,jci1,jci2,ici1,ici2,'surf:tg')
      call getmem2d(sfs%hfx,jci1,jci2,ici1,ici2,'surf:hfx')
      call getmem2d(sfs%qfx,jci1,jci2,ici1,ici2,'surf:qfx')
      call getmem2d(sfs%rainc,jci1,jci2,ici1,ici2,'surf:rainc')
      call getmem2d(sfs%rainnc,jci1,jci2,ici1,ici2,'surf:rainnc')
      if ( ipptls > 1 ) then
        call getmem2d(sfs%snownc,jci1,jci2,ici1,ici2,'surf:snownc')
      end if
      call getmem2d(sfs%tgbb,jci1,jci2,ici1,ici2,'surf:tgbb')
      call getmem2d(sfs%uvdrag,jci1,jci2,ici1,ici2,'surf:uvdrag')
      call getmem2d(sfs%zo,jci1,jci2,ici1,ici2,'surf:zo')
      call getmem2d(sfs%ram1,jci1,jci2,ici1,ici2,'surf:ram1')
      call getmem2d(sfs%rah1,jci1,jci2,ici1,ici2,'surf:rah1')
      call getmem2d(sfs%br,jci1,jci2,ici1,ici2,'surf:br')
      call getmem2d(sfs%q2m,jci1,jci2,ici1,ici2,'surf:q2m')
      call getmem2d(sfs%ustar,jci1,jci2,ici1,ici2,'surf:ustar')
      call getmem2d(sfs%w10m,jci1,jci2,ici1,ici2,'surf:w10m')
      call getmem2d(sfs%u10m,jci1,jci2,ici1,ici2,'surf:u10m')
      call getmem2d(sfs%v10m,jci1,jci2,ici1,ici2,'surf:v10m')
      if ( ibltyp == 4 ) then
        call getmem2d(sfs%uz0,jci1,jci2,ici1,ici2,'surf:uz0')
        call getmem2d(sfs%vz0,jci1,jci2,ici1,ici2,'surf:vz0')
        call getmem2d(sfs%thz0,jci1,jci2,ici1,ici2,'surf:thz0')
        call getmem2d(sfs%qz0,jci1,jci2,ici1,ici2,'surf:qz0')
      end if
    end subroutine allocate_surfstate

    subroutine allocate_slice(ax,a0)
      implicit none
      type(slice) , intent(out) :: ax
      type(reference_atmosphere) , intent(in) :: a0
      if ( idynamic == 3 ) then
        call getmem3d(ax%pf3d,jce1,jce2,ice1,ice2,1,kzp1,'slice:pf3d')
        call getmem3d(ax%zq,jce1,jce2,ice1,ice2,1,kzp1,'slice:zq')
        call getmem3d(ax%dzq,jce1,jce2,ice1,ice2,1,kz,'slice:dzq')
        call getmem3d(ax%rhb3d,jci1,jci2,ici1,ici2,1,kz,'slice:rhb3d')
        if ( icldmstrat == 1 ) then
          call getmem2d(ax%th700,jci1,jci2,ici1,ici2,'slice:th700')
        end if
        if ( ibltyp == 4 ) then
          call getmem3d(ax%tkepbl,jci1,jci2,ici1,ici2,1,kz,'slice:tkepbl')
        end if
      else
        call getmem3d(ax%ubx3d,jce1,jce2,ice1,ice2,1,kz,'slice:ubx3d')
        call getmem3d(ax%vbx3d,jce1,jce2,ice1,ice2,1,kz,'slice:vbx3d')
        call getmem3d(ax%pf3d,jce1,jce2,ice1,ice2,1,kzp1,'slice:pf3d')
        call getmem3d(ax%pb3d,jce1,jce2,ice1,ice2,1,kz,'slice:pb3d')
        call getmem3d(ax%rhob3d,jci1,jci2,ici1,ici2,1,kz,'slice:rhob3d')
        call getmem3d(ax%qsb3d,jce1,jce2,ice1,ice2,1,kz,'slice:qsb3d')
        call getmem3d(ax%rhb3d,jci1,jci2,ici1,ici2,1,kz,'slice:rhb3d')
        call getmem3d(ax%tv3d,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'slice:tv3d')
        if ( icldmstrat == 1 ) then
          call getmem2d(ax%th700,jci1,jci2,ici1,ici2,'slice:th700')
        end if
        call getmem3d(ax%ubd3d,jd1,jd2,id1,id2,1,kz,'slice:ubd3d')
        call getmem3d(ax%vbd3d,jd1,jd2,id1,id2,1,kz,'slice:vbd3d')
        call getmem3d(ax%tb3d,jx1,jx2,ix1,ix2,1,kz,'slice:tb3d')
        call getmem4d(ax%qxb3d,jx1,jx2,ix1,ix2,1,kz,1,nqx,'slice:qxb3d')
        if ( idynamic == 2 ) then
          call getmem3d(ax%ppb3d,jx1,jx2,ix1,ix2,1,kzp1,'slice:ppb3d')
          call getmem3d(ax%wb3d,jx1,jx2,ix1,ix2,1,kzp1,'slice:wb3d')
        end if
        if ( ichem == 1 ) then
          call getmem4d(ax%chib3d,jx1,jx2,ix1,ix2,1,kz,1,ntr,'slice:chib3d')
        end if
        if ( idynamic == 1 ) then
          call getmem3d(ax%zq,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'slice:zq')
          call getmem3d(ax%za,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'slice:za')
          call getmem3d(ax%dzq,jce1,jce2,ice1,ice2,1,kz,'slice:dzq')
          call getmem3d(ax%wb3d,jce1,jce2,ice1,ice2,1,kzp1,'slice:wb3d')
          call getmem3d(ax%wpx3d,jci1,jci2,ici1,ici2,1,kz,'slice:wpx3d')
        else
          call assignpnt(a0%z,atms%za)
          call assignpnt(a0%zf,atms%zq)
          call assignpnt(a0%dzf,atms%dzq)
        end if
        call getmem2d(ax%ps2d,jce1,jce2,ice1,ice2,'slice:ps2d')
        if ( ibltyp == 4 ) then
          call getmem3d(ax%tkepbl,jci1,jci2,ici1,ici2,1,kz,'slice:tkepbl')
        end if
      end if
      call getmem3d(ax%th3d,jce1,jce2,ice1,ice2,1,kz,'slice:th3d')
      call getmem2d(ax%rhox2d,jci1,jci2,ici1,ici2,'slice:rhox2d')
      call getmem2d(ax%tp2d,jci1,jci2,ici1,ici2,'slice:tp2d')
    end subroutine allocate_slice

    subroutine allocate_nhbh(nhbh)
      implicit none
      type(nhboundhelp) , intent(out) :: nhbh

      call getmem2d(nhbh%ps,jce1,jce2,ice1,ice2,'nhboundhelp:ps')
      if ( ichem == 1 .or. iclimaaer == 1 ) then
        call getmem3d(nhbh%tvirt,jce1,jce2,ice1,ice2,1,kz,'nhboundhelp:tvirt')
      end if
    end subroutine allocate_nhbh

    subroutine allocate_mod_atm_interface
      implicit none

      if ( idiffu == 1 ) then
        id1 = ide1gb
        id2 = ide2gb
        jd1 = jde1gb
        jd2 = jde2gb
        ix1 = ice1gb
        ix2 = ice2gb
        jx1 = jce1gb
        jx2 = jce2gb
      else if ( idiffu == 2 ) then
        id1 = ide1ga
        id2 = ide2ga
        jd1 = jde1ga
        jd2 = jde2ga
        ix1 = ice1ga
        ix2 = ice2ga
        jx1 = jce1ga
        jx2 = jce2ga
      else if ( idiffu == 3 ) then
        id1 = ide1gc
        id2 = ide2gc
        jd1 = jde1gc
        jd2 = jde2gc
        ix1 = ice1gc
        ix2 = ice2gc
        jx1 = jce1gc
        jx2 = jce2gc
      end if

      call allocate_domain(mddom)
      call allocate_domain_subgrid(mdsub)

      if ( idynamic == 2 ) then
        call allocate_reference_atmosphere(atm0)
        call allocate_nhbh(nhbh0)
        call allocate_nhbh(nhbh1)
      end if
      if ( idynamic == 3 ) then
        call allocate_atmosphere(mo_atm)
        call allocate_nhbh(nhbh0)
        call allocate_nhbh(nhbh1)
      else
        call allocate_atmstate_a(atm1)
        call allocate_atmstate_b(atm2)
        call allocate_atmstate_decoupled(atmx)
        call allocate_atmstate_c(atmc)
        call allocate_atmstate_tendency(aten)
        call allocate_mass_divergence(mdv)
      end if

      call allocate_slice(atms,atm0)
      call allocate_surfstate(sfs)

      ! FAB:
      !    complete for diag on water quantitiies idiag = 2, 3 etc
      if ( idiag > 0 ) then
        call allocate_tendiag(tdiag)
        call allocate_qendiag(qdiag)
      end if

      if ( idynamic == 1 ) then
        call getmem3d(dstor,jde1,jde2,ide1,ide2,1,nsplit,'storage:dstor')
        call getmem3d(hstor,jde1,jde2,ide1,ide2,1,nsplit,'storage:hstor')
      end if

      call getmem3d(omega,jci1,jci2,ici1,ici2,1,kz,'storage:omega')
      call getmem3d(qdot,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'storage:qdot')
      call getmem2d(ktrop,jci1,jci2,ici1,ici2,'storage:ktrop')
      call getmem2d(coszrs,jci1,jci2,ici1,ici2,'storage:coszrs')
      call getmem2d(pptc,jci1,jci2,ici1,ici2,'storage:pptc')
      call getmem2d(prca,jci1,jci2,ici1,ici2,'storage:prca')
      call getmem2d(icumbot,jci1,jci2,ici1,ici2,'storage:icumbot')
      call getmem2d(icumtop,jci1,jci2,ici1,ici2,'storage:icumtop')

      ! This needs to be saved in SAV file
      call getmem3d(fcc,jci1,jci2,ici1,ici2,1,kz,'storage:fcc')
      if ( ichem == 1 ) then
        call getmem3d(rembc,jci1,jci2,ici1,ici2,1,kz,'storage:rembc')
        call getmem3d(remrat,jci1,jci2,ici1,ici2,1,kz,'storage:remrat')
        call getmem3d(convpr,jci1,jci2,ici1,ici2,1,kz,'storage:convpr')
        call getmem2d(ssw2da,jci1,jci2,ici1,ici2,'storage:ssw2da')
        call getmem2d(sfracv2d,jci1,jci2,ici1,ici2,'storage:sfracv2d')
        call getmem2d(sfracb2d,jci1,jci2,ici1,ici2,'storage:sfracb2d')
        call getmem2d(sfracs2d,jci1,jci2,ici1,ici2,'storage:sfracs2d')
        call getmem2d(svegfrac2d,jci1,jci2,ici1,ici2,'storage:svegfrac2d')
        call getmem2d(sxlai2d,jci1,jci2,ici1,ici2,'storage:sxlai2d')
#ifdef CLM45
        call getmem3d(voc_em_clm,jci1,jci2, &
                                 ici1,ici2,1,ntr,'storage:voc_em_clm')
        call getmem3d(dustflx_clm,jci1,jci2, &
                                  ici1,ici2,1,4,'storage:dustflx_clm')
        call getmem3d(sw_vol,jci1,jci2, &
                             ici1,ici2,1,num_soil_layers,'storage:sw_vol')
        call getmem3d(tsoi,jci1,jci2, &
                           ici1,ici2,1,num_soil_layers,'storage:tsoi')
#endif
        call getmem3d(drydepflx,jci1,jci2,ici1,ici2,1,ntr,'storage:drydepflx')
        call getmem3d(wetdepflx,jci1,jci2,ici1,ici2,1,ntr,'storage:wetdepflx')
        if ( iindirect > 0 .and. iaerosol == 1 ) then
          call getmem3d(ccn,jci1,jci2,ici1,ici2,1,kz,'storage:ccn')
        end if
      end if
      if ( iocncpl == 1 .or. iwavcpl == 1 ) then
        call getmem2d(cplmsk,jci1,jci2,ici1,ici2,'storage:cplmsk')
        cplmsk(:,:) = 0
        call getmem3d(dailyrnf,jci1,jci2,ici1,ici2,1,2,'storage:dailyrnf')
      end if
      call getmem2d(pptnc,jci1,jci2,ici1,ici2,'storage:pptnc')
      call getmem2d(prnca,jci1,jci2,ici1,ici2,'storage:prnca')
      call getmem2d(ptrop,jci1,jci2,ici1,ici2,'storage:ptrop')
      call getmem3d(cldfra,jci1,jci2,ici1,ici2,1,kz,'storage:cldfra')
      call getmem3d(cldlwc,jci1,jci2,ici1,ici2,1,kz,'storage:cldlwc')
      call getmem3d(heatrt,jci1,jci2,ici1,ici2,1,kz,'storage:heatrt')
      call getmem2d(totcf,jci1,jci2,ici1,ici2,'storage:totcf')
      call getmem2d(flw,jci1,jci2,ici1,ici2,'storage:flw')
      call getmem2d(flwd,jci1,jci2,ici1,ici2,'storage:flwd')
      call getmem2d(fsw,jci1,jci2,ici1,ici2,'storage:fsw')
      call getmem2d(sabveg,jci1,jci2,ici1,ici2,'storage:sabveg')
      call getmem2d(solis,jci1,jci2,ici1,ici2,'storage:solis')
      call getmem2d(dsol,jci1,jci2,ici1,ici2,'storage:dsol')
      call getmem2d(solvs,jci1,jci2,ici1,ici2,'storage:solvs')
      call getmem2d(solvsd,jci1,jci2,ici1,ici2,'storage:solvsd')
      call getmem2d(solvl,jci1,jci2,ici1,ici2,'storage:solvl')
      call getmem2d(solvld,jci1,jci2,ici1,ici2,'storage:solvld')
      call getmem2d(albvl,jci1,jci2,ici1,ici2,'storage:albvl')
      call getmem2d(albvs,jci1,jci2,ici1,ici2,'storage:albvs')
      call getmem2d(aldirs,jci1,jci2,ici1,ici2,'storage:aldirs')
      call getmem2d(aldifs,jci1,jci2,ici1,ici2,'storage:aldifs')
      call getmem2d(aldirl,jci1,jci2,ici1,ici2,'storage:aldirl')
      call getmem2d(aldifl,jci1,jci2,ici1,ici2,'storage:aldifl')
      call getmem2d(emiss,jci1,jci2,ici1,ici2,'storage:emiss')
      call getmem2d(sinc,jci1,jci2,ici1,ici2,'storage:sinc')
      call getmem2d(sdelq,jci1,jci2,ici1,ici2,'storage:sdelq')
      call getmem2d(sdelt,jci1,jci2,ici1,ici2,'storage:sdelt')
      call getmem2d(zpbl,jci1,jci2,ici1,ici2,'storage:zpbl')
      call getmem2d(kpbl,jci1,jci2,ici1,ici2,'storage:kpbl')
      call getmem3d(q_detr,jci1,jci2,ici1,ici2,1,kz,'storage:q_detr')
      if ( any(icup == 5) ) then
        call getmem3d(rain_cc,jci1,jci2,ici1,ici2,1,kz,'storage:rain_cc')
      end if

      if ( ipptls == 2 ) then
        call getmem3d(rain_ls,jci1,jci2,ici1,ici2,1,kzp1,'storage:rain_ls')
      end if
      if ( idynamic == 2 ) then
        call getmem2d(dpsdxm,jce1,jce2,ice1,ice2,'storage:dpsdxm')
        call getmem2d(dpsdym,jce1,jce2,ice1,ice2,'storage:dpsdym')
      end if
      call getmem2d(crrate,jci1,jci2,ici1,ici2,'storage:crrate')
      call getmem2d(ncrrate,jci1,jci2,ici1,ici2,'storage:ncrrate')
    end subroutine allocate_mod_atm_interface

    subroutine export_data_from_atm(expfie)
      implicit none
      type(exp_data3d) , intent(inout) :: expfie
      integer(ik4) :: k , j , i

      call exchange(atm1%u,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm1%v,1,jde1,jde2,ide1,ide2,1,kz)

      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            expfie%u(j,i,k) = d_rfour*(atm1%u(j,i,k)+atm1%u(j+1,i,k) + &
                              atm1%u(j,i+1,k)+atm1%u(j+1,i+1,k)) /     &
                              sfs%psa(j,i)
            expfie%v(j,i,k) = d_rfour*(atm1%v(j,i,k)+atm1%v(j+1,i,k) + &
                              atm1%v(j,i+1,k)+atm1%v(j+1,i+1,k)) /     &
                              sfs%psa(j,i)
            if ( idynamic == 2 ) then
              expfie%w(j,i,k) = d_half*(atm1%w(j,i,k+1)+atm1%w(j,i,k)) / &
                                sfs%psa(j,i)
            else
              expfie%w(j,i,k) = (-1.0d0*omega(j,i,k)*d_1000) / &
                                (atm1%rho(j,i,k)*egrav)
            end if
            expfie%t(j,i,k) = atm1%t(j,i,k)/sfs%psa(j,i)
            !expfie%q(j,i,k) = atm1%qx(j,i,k,iqv)/sfs%psa(j,i)
            expfie%q(j,i,k) = atms%rhb3d(j,i,k)
          end do
        end do
      end do
    end subroutine export_data_from_atm

end module mod_atm_interface

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
