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

  implicit none (type, external)

  private

  logical, public, parameter :: cross = .false.
  logical, public, parameter :: dot = .true.

  type(domain), public :: mddom
  type(domain_subgrid), public :: mdsub
  type(atmstate_a), public :: atm1
  type(atmstate_b), public :: atm2
  type(atmstate_c), public :: atmc
  type(atmstate_tendency), public :: aten
  type(atmstate_decoupled), public :: atmx
  type(tendiag), public :: tdiag
  type(qendiag), public :: qdiag
  type(surfstate), public :: sfs
  type(slice), public :: atms
  type(v3dbound), public :: xtb, xqb, xub, xvb, xppb, xwwb, xpaib
  type(v3dbound), public :: xlb, xib
  type(v2dbound), public :: xpsb, xtsb
  type(bound_area), public :: ba_cr, ba_dt, ba_ut, ba_vt
  type(reference_atmosphere), public :: atm0
  type(mass_divergence), public :: mdv
  type(nhboundhelp), public :: nhbh0, nhbh1

  ! For idynamic 3

  type(atmosphere), public :: mo_atm

  public :: allocate_mod_atm_interface
  public :: allocate_v3dbound, allocate_v2dbound
  public :: setup_boundaries, setup_model_indexes
  public :: export_data_from_atm

  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: dstor => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: hstor => null( )

  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: qdot => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: omega => null( )

  ! Sun
  ! Cosine of zenithal solar angle
  real(rkx), public, pointer, contiguous, dimension(:,:) :: coszrs => null( )

  ! Cumulus
  integer(ik4), public, pointer, contiguous, dimension(:,:) :: icumbot => null( )
  integer(ik4), public, pointer, contiguous, dimension(:,:) :: icumtop => null( )
  integer(ik4), public, pointer, contiguous, dimension(:,:) :: ktrop => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: convpr => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:) :: pptc => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:) :: sptc => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:) :: prca => null( )

  ! Radiation
  real(rkx), pointer, contiguous, public, dimension(:,:) :: ptrop => null( )
  ! vegetation absorbed radiation (full solar spectrum)
  real(rkx), pointer, contiguous, public, dimension(:,:) :: sabveg => null( )
  ! Incident solar flux
  real(rkx), pointer, contiguous, public, dimension(:,:) :: dsol => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: solis => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: solvs => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: solvsd => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: solvl => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: solvld => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: totcf => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: flw => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: fsw => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: flwd => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: cldfra => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: cldlwc => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: heatrt => null( )

  ! Dynamic 2
  real(rkx), pointer, contiguous, public, dimension(:,:) :: dpsdxm => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: dpsdym => null( )
  real(rkx), public, dimension(-6:6,-6:6) :: tmask

  ! Surface
  ! Total Long wave albedo (0.7-5.0 micro-meter)
  real(rkx), pointer, contiguous, public, dimension(:,:) :: albvl => null( )
  ! Total Short wave albedo (0.2-0.7 micro-meter)
  real(rkx), pointer, contiguous, public, dimension(:,:) :: albvs => null( )
  ! 0.2-0.7 micro-meter srfc alb to direct radiation
  real(rkx), pointer, contiguous, public, dimension(:,:) :: aldirs => null( )
  ! 0.2-0.7 micro-meter srfc alb to diffuse radiation
  real(rkx), pointer, contiguous, public, dimension(:,:) :: aldifs => null( )
  ! 0.7-5.0 micro-meter srfc alb to direct radiation
  real(rkx), pointer, contiguous, public, dimension(:,:) :: aldirl => null( )
  ! 0.7-5.0 micro-meter srfc alb to diffuse radiation
  real(rkx), pointer, contiguous, public, dimension(:,:) :: aldifl => null( )
  ! Emissivity at surface
  real(rkx), pointer, contiguous, public, dimension(:,:) :: emiss => null( )
  ! Total solar incoming radiation
  real(rkx), pointer, contiguous,  public, dimension(:,:) :: sinc => null( )

  ! Precip
  real(rkx), pointer, contiguous, public, dimension(:,:) :: pptnc => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: prnca => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: crrate => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: ncrrate => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: fcc => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: remrat => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: rembc => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: ccn => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:,:) :: rain_ls => null( )

  ! PBL
  integer(ik4), public, pointer, contiguous, dimension(:,:) :: kpbl => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:) :: zpbl => null( )

  ! Cumulus
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: q_detr => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: rain_cc => null( )

  ! Surface for chemistry
  real(rkx), pointer, contiguous, public, dimension(:,:) :: sdelq => null( )
  real(rkx), pointer, contiguous, public, dimension(:,:) :: sdelt => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:) :: ssw2da => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:) :: sfracv2d => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:) :: sfracb2d => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:) :: sfracs2d => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:) :: svegfrac2d => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:) :: sxlai2d => null( )

#ifdef CLM45
  ! real(rkx), public, pointer, contiguous, dimension(:,:) :: ustar => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: voc_em_clm => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: dustflx_clm => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: ddepv_clm => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: sw_vol => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: tsoi => null( )
#endif

  !chemistry for surface
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: wetdepflx => null( )
  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: drydepflx => null( )

  ! Coupling
  integer(ik4), public, pointer, contiguous, dimension(:,:) :: cplmsk => null( )

  integer(ik4) :: ix1, ix2, jx1, jx2
  integer(ik4) :: id1, id2, jd1, jd2

#ifdef DEBUG
  !type(grid_nc_var4d), public, save :: nc_4d
  !type(grid_nc_var3d), public, save :: nc_3d
  !type(grid_nc_var2d), public, save :: nc_2d
  !type(grid_nc_var4d), public, save :: qqxp
#endif

  contains

    subroutine setup_model_indexes
      implicit none (type, external)
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
      write(ndebug,*) 'DOTPEXTJI1 : ', jde1, jde2, ide1, ide2
      write(ndebug,*) 'DOTPINTJI1 : ', jdi1, jdi2, idi1, idi2
      write(ndebug,*) 'DOTPINTJI2 : ', jdii1, jdii2, idii1, idii2
      write(ndebug,*) 'CRXPEXTJI1 : ', jce1, jce2, ice1, ice2
      write(ndebug,*) 'CRXPINTJI1 : ', jci1, jci2, ici1, ici2
      write(ndebug,*) 'CRXPINTJI2 : ', jcii1, jcii2, icii1, icii2
      write(ndebug,*) 'TOPBDY   : ', ma%has_bdytop
      write(ndebug,*) 'BTMBDY   : ', ma%has_bdybottom
      write(ndebug,*) 'RGTBDY   : ', ma%has_bdyright
      write(ndebug,*) 'LFTBDY   : ', ma%has_bdyleft
      flush(ndebug)
#endif
    end subroutine setup_model_indexes

    subroutine setup_boundaries(ldotx,ldoty,ba)
      implicit none (type, external)
      logical, intent(in) :: ldotx, ldoty
      type(bound_area), intent(inout) :: ba
      integer(ik4) :: icx, icy
      integer(ik4) :: igbb1, igbb2, igbt1, igbt2
      integer(ik4) :: jgbl1, jgbl2, jgbr1, jgbr2
      integer(ik4) :: i, j, i1, i2, j1, j2

      call getmem(ba%ibnd,jde1,jde2,ide1,ide2,'setup_boundaries:ibnd')
      call getmem(ba%bsouth,jde1,jde2,ide1,ide2,'setup_boundaries:bsouth')
      call getmem(ba%bnorth,jde1,jde2,ide1,ide2,'setup_boundaries:bnorth')
      call getmem(ba%beast,jde1,jde2,ide1,ide2,'setup_boundaries:beast')
      call getmem(ba%bwest,jde1,jde2,ide1,ide2,'setup_boundaries:bwest')

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
          do i = i1, i2
            if ( i >= igbb1 .and. i <= igbb2 ) then
              do j = j1, j2
                if ( j < jgbl1 .and. j > jgbr2 ) cycle
                ba%ibnd(j,i) = i-igbb1+2
                ba%bsouth(j,i) = .true.
                ba%ns = ba%ns+1
              end do
            end if
          end do
          ! North Boundary
          do i = i1, i2
            if ( i >= igbt1 .and. i <= igbt2 ) then
              do j = j1, j2
                if ( j < jgbl1 .and. j > jgbr2 ) cycle
                ba%ibnd(j,i) = igbt2-i+2
                ba%bnorth(j,i) = .true.
                ba%nn = ba%nn+1
              end do
            end if
          end do
        else
          ! Check for South boundary
          do i = i1, i2
            if ( i >= igbb1 .and. i <= igbb2 ) then
              do j = j1, j2
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
          do i = i1, i2
            if ( i >= igbt1 .and. i <= igbt2 ) then
              do j = j1, j2
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
          do i = i1, i2
            if ( i < igbb1 .or. i > igbt2 ) cycle
            do j = j1, j2
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
          do i = i1, i2
            if ( i < igbb1 .or. i > igbt2 ) cycle
            do j = j1, j2
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

      do i = ide2, ide1, -1
        do j = jde1, jde2
          if ( ba%bsouth(j,i) ) then
            write(ndebug,'(1a,i0.4)',advance='no') 'S', ba%ibnd(j,i)
          else if ( ba%bnorth(j,i) ) then
            write(ndebug,'(1a,i0.4)',advance='no') 'N', ba%ibnd(j,i)
          else if ( ba%bwest(j,i) ) then
            write(ndebug,'(1a,i0.4)',advance='no') 'W', ba%ibnd(j,i)
          else if ( ba%beast(j,i) ) then
            write(ndebug,'(1a,i0.4)',advance='no') 'E', ba%ibnd(j,i)
          else
            write(ndebug,'(1a,i0.4)',advance='no') 'X', 0
          end if
        end do
        write(ndebug,*) ' '
      end do
#endif
    end subroutine setup_boundaries

    subroutine allocate_v3dbound(xb,ke,ldot)
      implicit none (type, external)
      type(v3dbound), intent(inout) :: xb
      integer(ik4), intent(in) :: ke
      logical, intent(in) :: ldot
      if ( ldot ) then
        call getmem(xb%b0,jde1ga,jde2ga,ide1ga,ide2ga,1,ke,'v3dbound:b0')
        call getmem(xb%b1,jde1ga,jde2ga,ide1ga,ide2ga,1,ke,'v3dbound:b1')
        call getmem(xb%bt,jde1ga,jde2ga,ide1ga,ide2ga,1,ke,'v3dbound:bt')
      else
        call getmem(xb%b0,jce1ga,jce2ga,ice1ga,ice2ga,1,ke,'v3dbound:b0')
        call getmem(xb%b1,jce1ga,jce2ga,ice1ga,ice2ga,1,ke,'v3dbound:b1')
        call getmem(xb%bt,jce1ga,jce2ga,ice1ga,ice2ga,1,ke,'v3dbound:bt')
      end if
    end subroutine allocate_v3dbound

    subroutine allocate_v2dbound(xb,ldot)
      implicit none (type, external)
      type(v2dbound), intent(inout) :: xb
      logical, intent(in) :: ldot
      if ( ldot ) then
        call getmem(xb%b0,jde1ga,jde2ga,ide1ga,ide2ga,'v2dbound:b0')
        call getmem(xb%b1,jde1ga,jde2ga,ide1ga,ide2ga,'v2dbound:b1')
        call getmem(xb%bt,jde1ga,jde2ga,ide1ga,ide2ga,'v2dbound:bt')
      else
        call getmem(xb%b0,jce1ga,jce2ga,ice1ga,ice2ga,'v2dbound:b0')
        call getmem(xb%b1,jce1ga,jce2ga,ice1ga,ice2ga,'v2dbound:b1')
        call getmem(xb%bt,jce1ga,jce2ga,ice1ga,ice2ga,'v2dbound:bt')
      end if
    end subroutine allocate_v2dbound

    subroutine allocate_atmosphere(atm)
      implicit none (type, external)
      type(atmosphere), intent(inout) :: atm
      call getmem(atm%u,jde1gb,jde2gb,ice1ga,ice2ga,1,kz,'atmstate:u')
      call getmem(atm%v,jce1ga,jce2ga,ide1gb,ide2gb,1,kz,'atmstate:v')
#ifdef RCEMIP
      call getmem(atm%ux,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'atmstate:ux')
      call getmem(atm%vx,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'atmstate:vx')
#else
      call getmem(atm%ux,jce1gb,jce2gb,ice1ga,ice2ga,1,kz,'atmstate:ux')
      call getmem(atm%vx,jce1ga,jce2ga,ice1gb,ice2gb,1,kz,'atmstate:vx')
#endif
      call getmem(atm%t,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:t')
      call getmem(atm%tetav,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:tetav')
      call getmem(atm%w,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:w')
      call getmem(atm%pai,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:pai')
      call getmem(atm%qx,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,1,nqx,'atmstate:qx')
      call getmem(atm%zeta,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'atmstate:zeta')

      call getmem(atm%rho,jce1,jce2,ice1,ice2,1,kz,'atmstate:rho')
      call getmem(atm%pf,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:pf')
      call getmem(atm%p,jce1,jce2,ice1,ice2,1,kz,'atmstate:p')
      call getmem(atm%tvirt,jce1,jce2,ice1,ice2,1,kz,'atmstate:tvirt')
      call getmem(atm%zetaf,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:zetaf')
      call getmem(atm%dz,jce1,jce2,ice1,ice2,1,kz,'atmstate:dz')
      call getmem(atm%qs,jce1,jce2,ice1,ice2,1,kz,'atmstate:qs')

      call getmem(atm%tten,jci1,jci2,ici1,ici2,1,kz,'atmstate:tten')
      call getmem(atm%uten,jdi1,jdi2,ici1,ici2,1,kz,'atmstate:uten')
      call getmem(atm%vten,jci1,jci2,idi1,idi2,1,kz,'atmstate:vten')
      call getmem(atm%qxten,jci1,jci2,ici1,ici2,1,kz,1,nqx,'atmstate:qxten')
      if ( ibltyp == 2 ) then
        call getmem(atm%tke,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:tke')
        call getmem(atm%tketen,jci1,jci2,ici1,ici2,1,kzp1,'atmstate:tketen')
      end if
      if ( ichem == 1 ) then
        call getmem(atm%trac,jce1ga,jce2ga, &
                               ice1ga,ice2ga,1,kz,1,ntr,'atmstate:trac')
        call getmem(atm%chiten,jci1,jci2,   &
                                 ici1,ici2,1,kz,1,ntr,'atmstate:chiten')
      end if
      call getmem(atm%fmz,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:fmz')
      call getmem(atm%fmzf,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:fmzf')
    end subroutine allocate_atmosphere

    subroutine allocate_atmstate_a(atm)
      implicit none (type, external)
      type(atmstate_a), intent(inout) :: atm
      call getmem(atm%u,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:u')
      call getmem(atm%v,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:v')
      call getmem(atm%t,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:t')
      call getmem(atm%qx,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,1,nqx,'atmstate:qx')
      if ( ibltyp == 2 ) then
        call getmem(atm%tke,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem(atm%pr,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:pr')
        call getmem(atm%rho,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:rho')
        call getmem(atm%pp,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:pp')
        call getmem(atm%w,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'atmstate:w')
      else
        call getmem(atm%pr,jce1,jce2,ice1,ice2,1,kz,'atmstate:pr')
        call getmem(atm%rho,jce1,jce2,ice1,ice2,1,kz,'atmstate:rho')
      end if
      if ( ichem == 1 ) then
        call getmem(atm%chi,jce1ga,jce2ga,ice1ga,ice2ga, &
                              1,kz,1,ntr,'atmstate:chi')
      end if
    end subroutine allocate_atmstate_a

    subroutine allocate_atmstate_b(atm)
      implicit none (type, external)
      type(atmstate_b), intent(inout) :: atm

      call getmem(atm%u,jd1,jd2,id1,id2,1,kz,'atmstate:u')
      call getmem(atm%v,jd1,jd2,id1,id2,1,kz,'atmstate:v')
      call getmem(atm%t,jx1,jx2,ix1,ix2,1,kz,'atmstate:t')
      if ( isladvec == 1 ) then
        call getmem(atm%qx,jce1sl,jce2sl, &
                             ice1sl,ice2sl,1,kz,1,nqx,'atmstate:qx')
      else
        call getmem(atm%qx,jx1,jx2,ix1,ix2,1,kz,1,nqx,'atmstate:qx')
      end if
      if ( ibltyp == 2 ) then
        call getmem(atm%tke,jx1,jx2,ix1,ix2,1,kzp1,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem(atm%pp,jx1,jx2,ix1,ix2,1,kz,'atmstate:pp')
        call getmem(atm%w,jx1,jx2,ix1,ix2,1,kzp1,'atmstate:w')
      end if
      if ( ichem == 1 ) then
        if ( isladvec == 1 ) then
          call getmem(atm%chi,jce1sl,jce2sl,ice1sl,ice2sl, &
                                1,kz,1,ntr,'atmstate:chi')
        else
          call getmem(atm%chi,jx1,jx2,ix1,ix2,1,kz,1,ntr,'atmstate:chi')
        end if
      end if
      call getmem(atm%pr,jce1,jce2,ice1,ice2,1,kz,'atmstate:pr')
    end subroutine allocate_atmstate_b

    subroutine allocate_atmstate_c(atm)
      implicit none (type, external)
      type(atmstate_c), intent(inout) :: atm
      if ( ibltyp == 2 ) then
        call getmem(atm%tke,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem(atm%u,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:u')
        call getmem(atm%v,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:v')
        call getmem(atm%pp,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:pp')
        call getmem(atm%t,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:t')
        call getmem(atm%w,jce1,jce2,ice1,ice2,1,kzp1,'atmstate:w')
      else
        call getmem(atm%u,jdi1,jdi2,idi1,idi2,1,kz,'atmstate:u')
        call getmem(atm%v,jdi1,jdi2,idi1,idi2,1,kz,'atmstate:v')
        call getmem(atm%t,jci1,jci2,ici1,ici2,1,kz,'atmstate:t')
      end if
      call getmem(atm%qx,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,1,nqx,'atmstate:qx')
      if ( ichem == 1 ) then
        call getmem(atm%chi,jce1ga,jce2ga,ice1ga,ice2ga, &
                              1,kz,1,ntr,'atmstate:chi')
      end if
    end subroutine allocate_atmstate_c

    subroutine allocate_atmstate_decoupled(atm)
      implicit none (type, external)
      type(atmstate_decoupled), intent(inout) :: atm
      call getmem(atm%uc,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:uc')
      call getmem(atm%vc,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:vc')
      call getmem(atm%umc,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:umc')
      call getmem(atm%vmc,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:vmc')
      if ( isladvec == 1 ) then
        call getmem(atm%ud,jde1gb,jde2gb,ide1gb,ide2gb,1,kz,'atmstate:ud')
        call getmem(atm%vd,jde1gb,jde2gb,ide1gb,ide2gb,1,kz,'atmstate:vd')
        call getmem(atm%umd,jde1gb,jde2gb,ide1gb,ide2gb,1,kz,'atmstate:umd')
        call getmem(atm%vmd,jde1gb,jde2gb,ide1gb,ide2gb,1,kz,'atmstate:vmd')
      else
        call getmem(atm%ud,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:ud')
        call getmem(atm%vd,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:vd')
        call getmem(atm%umd,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:umd')
        call getmem(atm%vmd,jde1ga,jde2ga,ide1ga,ide2ga,1,kz,'atmstate:vmd')
      end if
      call getmem(atm%t,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:t')
      call getmem(atm%tv,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:tv')
      call getmem(atm%qx,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,1,nqx,'atmstate:qx')
      if ( idynamic == 2 ) then
        call getmem(atm%pr,jce1,jce2,ice1,ice2,1,kz,'atmstate:pr')
        call getmem(atm%pp,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'atmstate:pp')
        call getmem(atm%w,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'atmstate:w')
      end if
      if ( ichem == 1 ) then
        call getmem(atm%chi,jce1ga,jce2ga,ice1ga,ice2ga, &
                              1,kz,1,ntr,'atmstate:chi')
      end if
    end subroutine allocate_atmstate_decoupled

    subroutine allocate_atmstate_tendency(atm)
      implicit none (type, external)
      type(atmstate_tendency), intent(inout) :: atm
      call getmem(atm%u,jdi1,jdi2,idi1,idi2, &
                    1,kz,1,number_of_prognostic_components,'atmstate:u')
      call getmem(atm%v,jdi1,jdi2,idi1,idi2, &
                    1,kz,1,number_of_prognostic_components,'atmstate:v')
      call getmem(atm%t,jci1,jci2,ici1,ici2, &
                    1,kz,1,number_of_prognostic_components,'atmstate:t')
      call getmem(atm%qx,jci1,jci2,ici1,ici2, &
                    1,kz,1,nqx,1,number_of_prognostic_components,'atmstate:qx')
      if ( ibltyp == 2 ) then
        call getmem(atm%tke,jci1,jci2,ici1,ici2, &
                  1,kzp1,1,number_of_prognostic_components,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem(atm%pp,jci1,jci2,ici1,ici2, &
                      1,kz,1,number_of_prognostic_components,'atmstate:pp')
        call getmem(atm%w,jci1,jci2,ici1,ici2, &
                      1,kzp1,1,number_of_prognostic_components,'atmstate:w')
      end if
      if ( ichem == 1 ) then
        call getmem(atm%chi,jci1,jci2,ici1,ici2, &
                      1,kz,1,ntr,1,number_of_prognostic_components, &
                      'atmstate:chi')
      end if
    end subroutine allocate_atmstate_tendency

    subroutine allocate_reference_atmosphere(atm)
      implicit none (type, external)
      type(reference_atmosphere), intent(inout) :: atm
      call getmem(atm%ps,jce1ga,jce2ga,ice1ga,ice2ga,'reference:ps')
      call getmem(atm%psdot,jde1ga,jde2ga,ide1ga,ide2ga,'reference:psdot')
      call getmem(atm%t,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'reference:t')
      call getmem(atm%pr,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'reference:pr')
      call getmem(atm%rho,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'reference:rho')
      call getmem(atm%z,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'reference:z')
      call getmem(atm%zd,jdi1,jdi2,idi1,idi2,1,kz,'reference:zd')
      call getmem(atm%tf,jce1,jce2,ice1,ice2,1,kzp1,'reference:tf')
      call getmem(atm%pf,jce1,jce2,ice1,ice2,1,kzp1,'reference:pf')
      call getmem(atm%rhof,jce1,jce2,ice1,ice2,1,kzp1,'reference:rhof')
      call getmem(atm%zf,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'reference:zf')
      call getmem(atm%dzf,jce1,jce2,ice1,ice2,1,kzp1,'reference:dzf')
      call getmem(atm%dprddx,jdi1,jdi2,idi1,idi2,1,kz,'reference:dprddx')
      call getmem(atm%dprddy,jdi1,jdi2,idi1,idi2,1,kz,'reference:dprddy')
    end subroutine allocate_reference_atmosphere

    subroutine allocate_mass_divergence(div)
      implicit none (type, external)
      type(mass_divergence), intent(inout) :: div
      call getmem(div%cr,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'massdiv:cr')
    end subroutine allocate_mass_divergence

    subroutine allocate_tendiag(dia)
      implicit none (type, external)
      type(tendiag), intent(inout) :: dia
      call getmem(dia%adh,jci1,jci2,ici1,ici2,1,kz,'tendiag:adh')
      call getmem(dia%adv,jci1,jci2,ici1,ici2,1,kz,'tendiag:adv')
      call getmem(dia%tbl,jci1,jci2,ici1,ici2,1,kz,'tendiag:tbl')
      call getmem(dia%con,jci1,jci2,ici1,ici2,1,kz,'tendiag:con')
      call getmem(dia%bdy,jci1,jci2,ici1,ici2,1,kz,'tendiag:bdy')
      call getmem(dia%adi,jci1,jci2,ici1,ici2,1,kz,'tendiag:adi')
      call getmem(dia%dif,jci1,jci2,ici1,ici2,1,kz,'tendiag:dif')
      call getmem(dia%rad,jci1,jci2,ici1,ici2,1,kz,'tendiag:rad')
      call getmem(dia%lsc,jci1,jci2,ici1,ici2,1,kz,'tendiag:lsc')
    end subroutine allocate_tendiag

    subroutine allocate_qendiag(dia)
      implicit none (type, external)
      type(qendiag), intent(inout) :: dia
      call getmem(dia%adh,jci1,jci2,ici1,ici2,1,kz,'tendiag:adh')
      call getmem(dia%adv,jci1,jci2,ici1,ici2,1,kz,'tendiag:adv')
      call getmem(dia%tbl,jci1,jci2,ici1,ici2,1,kz,'tendiag:tbl')
      call getmem(dia%con,jci1,jci2,ici1,ici2,1,kz,'tendiag:con')
      call getmem(dia%bdy,jci1,jci2,ici1,ici2,1,kz,'tendiag:bdy')
      call getmem(dia%adi,jci1,jci2,ici1,ici2,1,kz,'tendiag:adi')
      call getmem(dia%dif,jci1,jci2,ici1,ici2,1,kz,'tendiag:dif')
      call getmem(dia%rad,jci1,jci2,ici1,ici2,1,kz,'tendiag:rad')
      call getmem(dia%lsc,jci1,jci2,ici1,ici2,1,kz,'tendiag:lsc')
      if ( ichem == 1 .and. iaerosol == 1 .and. iindirect == 2 ) then
        call getmem(dia%qcl,jci1,jci2,ici1,ici2,1,kz,'tendiag:qcl')
        call getmem(dia%qcr,jci1,jci2,ici1,ici2,1,kz,'tendiag:qcr')
        call getmem(dia%acr,jci1,jci2,ici1,ici2,1,kz,'tendiag:acr')
      end if
    end subroutine allocate_qendiag

    subroutine allocate_domain(dom)
      implicit none (type, external)
      type(domain), intent(inout) :: dom
      call getmem(dom%ht,jde1gb,jde2gb,ide1gb,ide2gb,'storage:ht')
      call getmem(dom%lndcat,jde1,jde2,ide1,ide2,'storage:lndcat')
      call getmem(dom%lndtex,jde1,jde2,ide1,ide2,'storage:lndtex')
      call getmem(dom%xlat,jde1ga,jde2ga,ide1ga,ide2ga,'storage:xlat')
      call getmem(dom%xlon,jde1ga,jde2ga,ide1ga,ide2ga,'storage:xlon')
      call getmem(dom%dlat,jde1,jde2,ide1,ide2,'storage:dlat')
      call getmem(dom%dlon,jde1,jde2,ide1,ide2,'storage:dlon')
      call getmem(dom%mask,jde1,jde2,ide1,ide2,'storage:mask')
      call getmem(dom%area,jde1,jde2,ide1,ide2,'storage:area')
      if ( idynamic == 3 ) then
        call getmem(dom%msfx,jde1,jde2,ide1,ide2,'storage:msfx')
        call getmem(dom%msfu,jde1ga,jde2ga,ide1,ide2,'storage:msfu')
        call getmem(dom%msfv,jde1,jde2,ide1ga,ide2ga,'storage:msfv')
        call getmem(dom%hx,jde1ga,jde2ga,ice1,ice2,'storage:hx')
        call getmem(dom%hy,jce1,jce2,ide1ga,ide2ga,'storage:hy')
        call getmem(dom%ulat,jde1,jde2,ide1,ide2,'storage:ulat')
        call getmem(dom%ulon,jde1,jde2,ide1,ide2,'storage:ulon')
        call getmem(dom%vlat,jde1,jde2,ide1,ide2,'storage:vlat')
        call getmem(dom%vlon,jde1,jde2,ide1,ide2,'storage:vlon')
        call getmem(dom%coriou,jde1,jde2,ice1,ice2,'storage:fu')
        call getmem(dom%coriov,jce1,jce2,ide1,ide2,'storage:fv')
      else
        call getmem(dom%msfx,jd1,jd2,id1,id2,'storage:msfx')
        call getmem(dom%msfd,jd1,jd2,id1,id2,'storage:msfd')
      end if
      call getmem(dom%coriol,jde1,jde2,ide1,ide2,'storage:f')
      call getmem(dom%snowam,jde1,jde2,ide1,ide2,'storage:snowam')
      call getmem(dom%smoist,jde1,jde2,ide1,ide2,'storage:smoist')
      call getmem(dom%rmoist,jde1,jde2,ide1,ide2, &
                    1,num_soil_layers,'storage:rmoist')
      call getmem(dom%rts,jde1,jde2,ide1,ide2, &
                    1,num_soil_layers,'storage:rts')
      call getmem(dom%ldmsk,jci1,jci2,ici1,ici2,'storage:ldmsk')
      call getmem(dom%iveg,jci1,jci2,ici1,ici2,'storage:iveg')
      call getmem(dom%itex,jci1,jci2,ici1,ici2,'storage:itex')
      call getmem(dom%xmsf,jdi1,jdi2,idi1,idi2,'storage:xmsf')
      call getmem(dom%dmsf,jdi1,jdi2,idi1,idi2,'storage:dmsf')
      if ( lakemod == 1 ) then
        call getmem(dom%dhlake,jde1,jde2,ide1,ide2,'storage:dhlake')
      end if
      if ( idynamic == 2 ) then
        call getmem(dom%ef,jdi1ga,jdi2ga,idi1ga,idi2ga,'storage:ef')
        call getmem(dom%ddx,jdi1ga,jdi2ga,idi1ga,idi2ga,'storage:ddx')
        call getmem(dom%ddy,jdi1ga,jdi2ga,idi1ga,idi2ga,'storage:ddy')
        call getmem(dom%ex,jci1,jci2,ici1,ici2,'storage:ex')
        call getmem(dom%crx,jci1,jci2,ici1,ici2,'storage:crx')
        call getmem(dom%cry,jci1,jci2,ici1,ici2,'storage:cry')
        call getmem(dom%dmdy,jdi1,jdi2,idi1,idi2,'storage:dmdy')
        call getmem(dom%dmdx,jdi1,jdi2,idi1,idi2,'storage:dmdx')
      end if
    end subroutine allocate_domain

    subroutine allocate_domain_subgrid(sub)
      implicit none (type, external)
      type(domain_subgrid), intent(inout) :: sub
      call getmem(sub%ht,1,nnsg,jde1,jde2,ide1,ide2,'storage:ht')
      call getmem(sub%lndcat,1,nnsg,jde1,jde2,ide1,ide2,'storage:lndcat')
      call getmem(sub%lndtex,1,nnsg,jde1,jde2,ide1,ide2,'storage:lndtex')
      call getmem(sub%xlat,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlat')
      call getmem(sub%xlon,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlon')
      call getmem(sub%mask,1,nnsg,jde1,jde2,ide1,ide2,'storage:mask')
      call getmem(sub%area,1,nnsg,jde1,jde2,ide1,ide2,'storage:area')
      call getmem(sub%ldmsk,1,nnsg,jci1,jci2,ici1,ici2,'storage:ldmsk')
      call getmem(sub%iveg,1,nnsg,jci1,jci2,ici1,ici2,'storage:iveg')
      call getmem(sub%itex,1,nnsg,jci1,jci2,ici1,ici2,'storage:itex')
      if ( lakemod == 1 ) then
        call getmem(sub%dhlake,1,nnsg,jde1,jde2,ide1,ide2,'storage:dhlake')
      end if
    end subroutine allocate_domain_subgrid

    subroutine allocate_surfstate(sfs)
      implicit none (type, external)
      type(surfstate), intent(inout) :: sfs
      call getmem(sfs%psa,jce1ga,jce2ga,ice1ga,ice2ga,'surf:psa')
      if ( idynamic == 3 ) then
        call assignpnt(sfs%psa,sfs%psb)
        call assignpnt(sfs%psa,sfs%psc)
      else
        call getmem(sfs%psdota,jde1ga,jde2ga,ide1ga,ide2ga,'surf:psdota')
        call getmem(sfs%psb,jx1,jx2,ix1,ix2,'surf:psb')
        call getmem(sfs%psdotb,jd1,jd2,id1,id2,'surf:psdotb')
        call getmem(sfs%psc,jce1,jce2,ice1,ice2,'surf:psc')
      end if
      call getmem(sfs%tg,jci1,jci2,ici1,ici2,'surf:tg')
      call getmem(sfs%hfx,jci1,jci2,ici1,ici2,'surf:hfx')
      call getmem(sfs%qfx,jci1,jci2,ici1,ici2,'surf:qfx')
      call getmem(sfs%rainc,jci1,jci2,ici1,ici2,'surf:rainc')
      call getmem(sfs%rainnc,jci1,jci2,ici1,ici2,'surf:rainnc')
      if ( ipptls > 1 ) then
        call getmem(sfs%snownc,jci1,jci2,ici1,ici2,'surf:snownc')
      end if
      if ( ipptls > 3 ) then
        call getmem(sfs%grplnc,jci1,jci2,ici1,ici2,'surf:grplnc')
        call getmem(sfs%hailnc,jci1,jci2,ici1,ici2,'surf:hailnc')
      end if
      call getmem(sfs%tgbb,jci1,jci2,ici1,ici2,'surf:tgbb')
      call getmem(sfs%uvdrag,jci1,jci2,ici1,ici2,'surf:uvdrag')
      call getmem(sfs%zo,jci1,jci2,ici1,ici2,'surf:zo')
      call getmem(sfs%ram1,jci1,jci2,ici1,ici2,'surf:ram1')
      call getmem(sfs%rah1,jci1,jci2,ici1,ici2,'surf:rah1')
      if ( iocncpl == 1 .or. iwavcpl == 1 ) then
        call getmem(sfs%dtrnof,jci1,jci2,ici1,ici2,'surf:dtrnof')
      end if
      call getmem(sfs%br,jci1,jci2,ici1,ici2,'surf:br')
      call getmem(sfs%q2m,jci1,jci2,ici1,ici2,'surf:q2m')
      call getmem(sfs%ustar,jci1,jci2,ici1,ici2,'surf:ustar')
      call getmem(sfs%w10m,jci1,jci2,ici1,ici2,'surf:w10m')
      call getmem(sfs%u10m,jci1,jci2,ici1,ici2,'surf:u10m')
      call getmem(sfs%v10m,jci1,jci2,ici1,ici2,'surf:v10m')
      if ( ibltyp == 4 ) then
        call getmem(sfs%uz0,jci1,jci2,ici1,ici2,'surf:uz0')
        call getmem(sfs%vz0,jci1,jci2,ici1,ici2,'surf:vz0')
        call getmem(sfs%thz0,jci1,jci2,ici1,ici2,'surf:thz0')
        call getmem(sfs%qz0,jci1,jci2,ici1,ici2,'surf:qz0')
      end if
    end subroutine allocate_surfstate

    subroutine allocate_slice(ax,a0)
      implicit none (type, external)
      type(slice), intent(inout) :: ax
      type(reference_atmosphere), intent(in) :: a0
      if ( idynamic == 3 ) then
        call getmem(ax%pf3d,jce1,jce2,ice1,ice2,1,kzp1,'slice:pf3d')
        call getmem(ax%zq,jce1,jce2,ice1,ice2,1,kzp1,'slice:zq')
        call getmem(ax%dzq,jce1,jce2,ice1,ice2,1,kz,'slice:dzq')
        call getmem(ax%rhb3d,jci1,jci2,ici1,ici2,1,kz,'slice:rhb3d')
        if ( icldmstrat == 1 ) then
          call getmem(ax%th700,jci1,jci2,ici1,ici2,'slice:th700')
        end if
        if ( ibltyp == 4 .or. ibltyp == 5 ) then
          call getmem(ax%tkepbl,jci1,jci2,ici1,ici2,1,kz,'slice:tkepbl')
        end if
      else
        call getmem(ax%ubx3d,jce1,jce2,ice1,ice2,1,kz,'slice:ubx3d')
        call getmem(ax%vbx3d,jce1,jce2,ice1,ice2,1,kz,'slice:vbx3d')
        call getmem(ax%pf3d,jce1,jce2,ice1,ice2,1,kzp1,'slice:pf3d')
        call getmem(ax%pb3d,jce1,jce2,ice1,ice2,1,kz,'slice:pb3d')
        call getmem(ax%rhob3d,jci1,jci2,ici1,ici2,1,kz,'slice:rhob3d')
        call getmem(ax%qsb3d,jce1,jce2,ice1,ice2,1,kz,'slice:qsb3d')
        call getmem(ax%rhb3d,jci1,jci2,ici1,ici2,1,kz,'slice:rhb3d')
        call getmem(ax%tv3d,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'slice:tv3d')
        if ( icldmstrat == 1 ) then
          call getmem(ax%th700,jci1,jci2,ici1,ici2,'slice:th700')
        end if
        call getmem(ax%ubd3d,jd1,jd2,id1,id2,1,kz,'slice:ubd3d')
        call getmem(ax%vbd3d,jd1,jd2,id1,id2,1,kz,'slice:vbd3d')
        call getmem(ax%tb3d,jx1,jx2,ix1,ix2,1,kz,'slice:tb3d')
        call getmem(ax%qxb3d,jx1,jx2,ix1,ix2,1,kz,1,nqx,'slice:qxb3d')
        if ( idynamic == 2 ) then
          call getmem(ax%ppb3d,jx1,jx2,ix1,ix2,1,kzp1,'slice:ppb3d')
          call getmem(ax%wb3d,jx1,jx2,ix1,ix2,1,kzp1,'slice:wb3d')
        end if
        if ( ichem == 1 ) then
          call getmem(ax%chib3d,jx1,jx2,ix1,ix2,1,kz,1,ntr,'slice:chib3d')
        end if
        if ( idynamic == 1 ) then
          call getmem(ax%zq,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'slice:zq')
          call getmem(ax%za,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'slice:za')
          call getmem(ax%dzq,jce1,jce2,ice1,ice2,1,kz,'slice:dzq')
          call getmem(ax%wb3d,jce1,jce2,ice1,ice2,1,kzp1,'slice:wb3d')
          call getmem(ax%wpx3d,jci1,jci2,ici1,ici2,1,kz,'slice:wpx3d')
        else
          call assignpnt(a0%z,atms%za)
          call assignpnt(a0%zf,atms%zq)
          call assignpnt(a0%dzf,atms%dzq)
        end if
        call getmem(ax%ps2d,jce1,jce2,ice1,ice2,'slice:ps2d')
        if ( ibltyp == 4 .or. ibltyp == 5 ) then
          call getmem(ax%tkepbl,jci1,jci2,ici1,ici2,1,kz,'slice:tkepbl')
        end if
      end if
      call getmem(ax%th3d,jce1,jce2,ice1,ice2,1,kz,'slice:th3d')
      call getmem(ax%rhox2d,jci1,jci2,ici1,ici2,'slice:rhox2d')
      call getmem(ax%tp2d,jci1,jci2,ici1,ici2,'slice:tp2d')
    end subroutine allocate_slice

    subroutine allocate_nhbh(nhbh)
      implicit none (type, external)
      type(nhboundhelp), intent(inout) :: nhbh

      call getmem(nhbh%ps,jce1,jce2,ice1,ice2,'nhboundhelp:ps')
      if ( ichem == 1 .or. iclimaaer == 1 ) then
        call getmem(nhbh%tvirt,jce1,jce2,ice1,ice2,1,kz,'nhboundhelp:tvirt')
      end if
    end subroutine allocate_nhbh

    subroutine allocate_mod_atm_interface
      implicit none (type, external)

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
        call getmem(dstor,jde1,jde2,ide1,ide2,1,nsplit,'storage:dstor')
        call getmem(hstor,jde1,jde2,ide1,ide2,1,nsplit,'storage:hstor')
      end if

      call getmem(omega,jci1,jci2,ici1,ici2,1,kz,'storage:omega')
      call getmem(qdot,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'storage:qdot')
      call getmem(ktrop,jci1,jci2,ici1,ici2,'storage:ktrop')
      call getmem(coszrs,jci1,jci2,ici1,ici2,'storage:coszrs')
      call getmem(prca,jci1,jci2,ici1,ici2,'storage:prca')
      call getmem(pptc,jci1,jci2,ici1,ici2,'storage:pptc')
      if ( ipptls > 1 .and. any(icup == 5) ) then
        call getmem(sptc,jci1,jci2,ici1,ici2,'storage:sptc')
      end if
      call getmem(icumbot,jci1,jci2,ici1,ici2,'storage:icumbot')
      call getmem(icumtop,jci1,jci2,ici1,ici2,'storage:icumtop')

      ! This needs to be saved in SAV file
      call getmem(fcc,jci1,jci2,ici1,ici2,1,kz,'storage:fcc')
      if ( ichem == 1 ) then
        call getmem(rembc,jci1,jci2,ici1,ici2,1,kz,'storage:rembc')
        call getmem(remrat,jci1,jci2,ici1,ici2,1,kz,'storage:remrat')
        call getmem(convpr,jci1,jci2,ici1,ici2,1,kz,'storage:convpr')
        call getmem(ssw2da,jci1,jci2,ici1,ici2,'storage:ssw2da')
        call getmem(sfracv2d,jci1,jci2,ici1,ici2,'storage:sfracv2d')
        call getmem(sfracb2d,jci1,jci2,ici1,ici2,'storage:sfracb2d')
        call getmem(sfracs2d,jci1,jci2,ici1,ici2,'storage:sfracs2d')
        call getmem(svegfrac2d,jci1,jci2,ici1,ici2,'storage:svegfrac2d')
        call getmem(sxlai2d,jci1,jci2,ici1,ici2,'storage:sxlai2d')
#ifdef CLM45
        call getmem(voc_em_clm,jci1,jci2, &
                                 ici1,ici2,1,ntr,'storage:voc_em_clm')
        call getmem(dustflx_clm,jci1,jci2, &
                                  ici1,ici2,1,4,'storage:dustflx_clm')
        call getmem(ddepv_clm,jci1,jci2, &
                                  ici1,ici2,1,ntr,'storage:ddepv_clm')
        call getmem(sw_vol,jci1,jci2, &
                             ici1,ici2,1,num_soil_layers,'storage:sw_vol')
        call getmem(tsoi,jci1,jci2, &
                           ici1,ici2,1,num_soil_layers,'storage:tsoi')
#endif
        call getmem(drydepflx,jci1,jci2,ici1,ici2,1,ntr,'storage:drydepflx')
        call getmem(wetdepflx,jci1,jci2,ici1,ici2,1,ntr,'storage:wetdepflx')
        if ( iindirect > 0 .and. iaerosol == 1 ) then
          call getmem(ccn,jci1,jci2,ici1,ici2,1,kz,'storage:ccn')
        end if
      end if
      if ( iocncpl == 1 .or. iwavcpl == 1 ) then
        call getmem(cplmsk,jci1,jci2,ici1,ici2,'storage:cplmsk')
        cplmsk(:,:) = 0
      end if
      call getmem(pptnc,jci1,jci2,ici1,ici2,'storage:pptnc')
      call getmem(prnca,jci1,jci2,ici1,ici2,'storage:prnca')
      call getmem(ptrop,jci1,jci2,ici1,ici2,'storage:ptrop')
      call getmem(cldfra,jci1,jci2,ici1,ici2,1,kz,'storage:cldfra')
      call getmem(cldlwc,jci1,jci2,ici1,ici2,1,kz,'storage:cldlwc')
      call getmem(heatrt,jci1,jci2,ici1,ici2,1,kz,'storage:heatrt')
      call getmem(totcf,jci1,jci2,ici1,ici2,'storage:totcf')
      call getmem(flw,jci1,jci2,ici1,ici2,'storage:flw')
      call getmem(flwd,jci1,jci2,ici1,ici2,'storage:flwd')
      call getmem(fsw,jci1,jci2,ici1,ici2,'storage:fsw')
      call getmem(sabveg,jci1,jci2,ici1,ici2,'storage:sabveg')
      call getmem(solis,jci1,jci2,ici1,ici2,'storage:solis')
      call getmem(dsol,jci1,jci2,ici1,ici2,'storage:dsol')
      call getmem(solvs,jci1,jci2,ici1,ici2,'storage:solvs')
      call getmem(solvsd,jci1,jci2,ici1,ici2,'storage:solvsd')
      call getmem(solvl,jci1,jci2,ici1,ici2,'storage:solvl')
      call getmem(solvld,jci1,jci2,ici1,ici2,'storage:solvld')
      call getmem(albvl,jci1,jci2,ici1,ici2,'storage:albvl')
      call getmem(albvs,jci1,jci2,ici1,ici2,'storage:albvs')
      call getmem(aldirs,jci1,jci2,ici1,ici2,'storage:aldirs')
      call getmem(aldifs,jci1,jci2,ici1,ici2,'storage:aldifs')
      call getmem(aldirl,jci1,jci2,ici1,ici2,'storage:aldirl')
      call getmem(aldifl,jci1,jci2,ici1,ici2,'storage:aldifl')
      call getmem(emiss,jci1,jci2,ici1,ici2,'storage:emiss')
      call getmem(sinc,jci1,jci2,ici1,ici2,'storage:sinc')
      call getmem(sdelq,jci1,jci2,ici1,ici2,'storage:sdelq')
      call getmem(sdelt,jci1,jci2,ici1,ici2,'storage:sdelt')
      call getmem(zpbl,jci1,jci2,ici1,ici2,'storage:zpbl')
      call getmem(kpbl,jci1,jci2,ici1,ici2,'storage:kpbl')
      call getmem(q_detr,jci1,jci2,ici1,ici2,1,kz,'storage:q_detr')
      if ( any(icup == 5) ) then
        call getmem(rain_cc,jci1,jci2,ici1,ici2,1,kz,'storage:rain_cc')
      end if

      if ( ipptls == 2 ) then
        call getmem(rain_ls,jci1,jci2,ici1,ici2,1,kz,'storage:rain_ls')
      end if
      if ( idynamic == 2 ) then
        call getmem(dpsdxm,jce1,jce2,ice1,ice2,'storage:dpsdxm')
        call getmem(dpsdym,jce1,jce2,ice1,ice2,'storage:dpsdym')
      end if
      call getmem(crrate,jci1,jci2,ici1,ici2,'storage:crrate')
      call getmem(ncrrate,jci1,jci2,ici1,ici2,'storage:ncrrate')
    end subroutine allocate_mod_atm_interface

    subroutine export_data_from_atm(expfie)
      implicit none (type, external)
      type(exp_data3d), intent(inout) :: expfie
      integer(ik4) :: k, j, i

      if ( idynamic == 3 ) then
        call exchange(mo_atm%u,1,jde1,jde2,ice1,ice2,1,kz)
        call exchange(mo_atm%v,1,jce1,jce2,ide1,ide2,1,kz)
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          expfie%u(j,i,k) = d_half*(mo_atm%u(j,i,k)+mo_atm%u(j+1,i,k))
          expfie%v(j,i,k) = d_half*(mo_atm%v(j,i,k)+mo_atm%u(j,i+1,k))
          expfie%w(j,i,k) = d_half*(mo_atm%w(j,i,k)+mo_atm%w(j,i,k+1))
          expfie%t(j,i,k) = mo_atm%t(j,i,k)
          expfie%q(j,i,k) = mo_atm%rho(j,i,k)
        end do
      else
        call exchange(atm1%u,1,jde1,jde2,ide1,ide2,1,kz)
        call exchange(atm1%v,1,jde1,jde2,ide1,ide2,1,kz)
        do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
          expfie%u(j,i,k) = d_rfour*(atm1%u(j,i,k)+atm1%u(j+1,i,k) + &
                            atm1%u(j,i+1,k)+atm1%u(j+1,i+1,k)) /     &
                            sfs%psa(j,i)
          expfie%v(j,i,k) = d_rfour*(atm1%v(j,i,k)+atm1%v(j+1,i,k) + &
                            atm1%v(j,i+1,k)+atm1%v(j+1,i+1,k)) /     &
                            sfs%psa(j,i)
          expfie%t(j,i,k) = atm1%t(j,i,k)/sfs%psa(j,i)
          !expfie%q(j,i,k) = atm1%qx(j,i,k,iqv)/sfs%psa(j,i)
          expfie%q(j,i,k) = atms%rhb3d(j,i,k)
        end do
        if ( idynamic == 2 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            expfie%w(j,i,k) = d_half*(atm1%w(j,i,k+1)+atm1%w(j,i,k)) / &
                              sfs%psa(j,i)
          end do
        else
          do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
            expfie%w(j,i,k) = (-1.0d0*omega(j,i,k)*d_1000) / &
                              (atm1%rho(j,i,k)*egrav)
          end do
        end if
      end if
    end subroutine export_data_from_atm

end module mod_atm_interface

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
