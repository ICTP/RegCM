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
  use mod_constants , only : d_rfour
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
  type(atmstate_tendency) , public :: aten , uwten
  type(atmstate_decoupled) , public :: atmx
  type(tendiag) , public :: tdiag
  type(qendiag) , public :: qdiag
  type(surfstate) , public :: sfs
  type(slice) , public :: atms
  type(diffx) , public :: adf
  type(v3dbound) , public :: xtb , xqb , xub , xvb
  type(v2dbound) , public :: xpsb
  type(bound_area) , public :: ba_cr , ba_dt
  type(reference_atmosphere) , public :: atm0

  public :: allocate_mod_atm_interface
  public :: allocate_v3dbound , allocate_v2dbound
  public :: setup_boundaries , setup_model_indexes

  real(rk8) , public , pointer , dimension(:,:) :: hgfact
  real(rk8) , public , pointer , dimension(:,:,:) :: dstor
  real(rk8) , public , pointer , dimension(:,:,:) :: hstor
  real(rk8) , public , pointer , dimension(:,:) :: ts0 , ts1

  real(rk8) , public , pointer , dimension(:,:,:) :: qdot , omega

  ! Sun
  ! Cosine of zenithal solar angle
  real(rk8) , public , pointer , dimension(:,:) :: coszrs

  ! Cumulus
  integer(ik4) , pointer , public , dimension(:,:) :: icumbot
  integer(ik4) , pointer , public , dimension(:,:) :: icumtop
  integer(ik4) , pointer , public , dimension(:,:) :: ktrop
  real(rk8) , pointer , public , dimension(:,:,:) :: convpr
  real(rk8) , pointer , public , dimension(:,:) :: pptc
  real(rk8) , pointer , public , dimension(:,:) :: prca

  ! Radiation
  real(rk8) , pointer , public , dimension(:,:) :: ptrop
  ! vegetation absorbed radiation (full solar spectrum)
  real(rk8) , pointer , public , dimension(:,:) :: sabveg
  ! Incident solar flux
  real(rk8) , pointer , public , dimension(:,:) :: solis
  real(rk8) , pointer , public , dimension(:,:) :: solvs
  real(rk8) , pointer , public , dimension(:,:) :: solvsd
  real(rk8) , pointer , public , dimension(:,:) :: solvl
  real(rk8) , pointer , public , dimension(:,:) :: solvld
  real(rk8) , pointer , public , dimension(:,:) :: flw
  real(rk8) , pointer , public , dimension(:,:) :: fsw
  real(rk8) , pointer , public , dimension(:,:) :: flwd
  real(rk8) , pointer , public , dimension(:,:,:) :: cldfra
  real(rk8) , pointer , public , dimension(:,:,:) :: cldlwc
  real(rk8) , pointer , public , dimension(:,:,:) :: heatrt

  ! Dynamic 2
  real(rk8) , pointer , public , dimension(:,:) :: dpsdxm , dpsdym

  ! Surface
  ! Total Long wave albedo (0.7-5.0 micro-meter)
  real(rk8) , pointer , public , dimension(:,:) :: albvl
  ! Total Short wave albedo (0.2-0.7 micro-meter)
  real(rk8) , pointer , public , dimension(:,:) :: albvs
  ! 0.2-0.7 micro-meter srfc alb to direct radiation
  real(rk8) , pointer , public , dimension(:,:) :: aldirs
  ! 0.2-0.7 micro-meter srfc alb to diffuse radiation
  real(rk8) , pointer , public , dimension(:,:) :: aldifs
  ! 0.7-5.0 micro-meter srfc alb to direct radiation
  real(rk8) , pointer , public , dimension(:,:) :: aldirl
  ! 0.7-5.0 micro-meter srfc alb to diffuse radiation
  real(rk8) , pointer , public , dimension(:,:) :: aldifl
  ! Emissivity at surface
  real(rk8) , pointer , public , dimension(:,:) :: emiss
  ! Total solar incoming radiation
  real(rk8) , pointer , public , dimension(:,:) :: sinc

  ! Precip
  real(rk8) , pointer , public , dimension(:,:) :: pptnc
  real(rk8) , pointer , public , dimension(:,:) :: prnca

  ! PBL
  integer(ik4) , public , pointer , dimension(:,:) :: kpbl
  real(rk8) , public , pointer , dimension(:,:) :: zpbl

  integer(ik4) , public , parameter :: zero_exchange_point = 0
  integer(ik4) , public , parameter :: one_exchange_point = 1
  integer(ik4) , public , parameter :: two_exchange_point = 2
  integer(ik4) , public , parameter :: four_exchange_point = 4

  ! Surface for chemistry
  real(rk8) , pointer , public , dimension(:,:) :: sdelq
  real(rk8) , pointer , public , dimension(:,:) :: sdelt
  real(rk8) , public , pointer , dimension(:,:) :: ssw2da
  real(rk8) , public , pointer , dimension(:,:) :: sfracv2d
  real(rk8) , public , pointer , dimension(:,:) :: sfracb2d
  real(rk8) , public , pointer , dimension(:,:) :: sfracs2d
  real(rk8) , public , pointer , dimension(:,:) :: svegfrac2d
  real(rk8) , public , pointer , dimension(:,:) :: sxlai2d

  ! Coupling
  real(rk8) , public , pointer , dimension(:,:,:) :: dailyrnf
  integer(ik4) , public , pointer , dimension(:,:) :: cplmsk

#ifdef DEBUG
  !type(grid_nc_var4d) , public , save :: nc_4d
  !type(grid_nc_var3d) , public , save :: nc_3d
  !type(grid_nc_var2d) , public , save :: nc_2d
  !type(grid_nc_var4d) , public , save :: qqxp
#endif

  contains
!
    subroutine setup_model_indexes
      implicit none
      ma%jbl1 = 1
      ma%jbl2 = 2
      ma%jbl4 = 4
      ma%jbr1 = 1
      ma%jbr2 = 2
      ma%jbr4 = 4
      ma%ibt1 = 1
      ma%ibt2 = 2
      ma%ibt4 = 4
      ma%ibb1 = 1
      ma%ibb2 = 2
      ma%ibb4 = 4
      if ( ma%has_bdyleft ) then
        ma%jbl1 = 0
        ma%jbl2 = 0
        ma%jbl4 = 0
      end if
      if ( ma%has_bdyright ) then
        ma%jbr1 = 0
        ma%jbr2 = 0
        ma%jbr4 = 0
      end if
      if ( ma%has_bdytop ) then
        ma%ibt1 = 0
        ma%ibt2 = 0
        ma%ibt4 = 0
      end if
      if ( ma%has_bdybottom ) then
        ma%ibb1 = 0
        ma%ibb2 = 0
        ma%ibb4 = 0
      end if
      jde1  = 1
      jdi1  = 1
      jdii1 = 1
      jde2  = jxp
      jdi2  = jxp
      jdii2 = jxp
      ide1  = 1
      idi1  = 1
      idii1 = 1
      ide2  = iyp
      idi2  = iyp
      idii2 = iyp
      if ( ma%has_bdyleft ) then
        jdi1 = 2
        jdii1 = 3
      end if
      if ( ma%has_bdyright ) then
        jdi2 = jxp-1
        jdii2 = jxp-2
      end if
      if ( ma%has_bdybottom ) then
        idi1 = 2
        idii1 = 3
      end if
      if ( ma%has_bdytop ) then
        idi2 = iyp-1
        idii2 = iyp-2
      end if
      jce1  = 1
      jci1  = 1
      jcii1 = 1
      jce2  = jxp
      jci2  = jxp
      jcii2 = jxp
      ice1  = 1
      ici1  = 1
      icii1 = 1
      ice2  = iyp
      ici2  = iyp
      icii2 = iyp
      if ( ma%has_bdyleft ) then
        jci1 = 2
        jcii1 = 3
      end if
      if ( ma%has_bdyright ) then
        jce2 = jxp-1
        jci2 = jxp-2
        jcii2 = jxp-3
      end if
      if ( ma%has_bdybottom ) then
        ici1 = 2
        icii1 = 3
      end if
      if ( ma%has_bdytop ) then
        ice2 = iyp-1
        ici2 = iyp-2
        icii2 = iyp-3
      end if
#ifdef DEBUG
      write(ndebug+myid,*) 'TOPLEFT     = ', ma%topleft
      write(ndebug+myid,*) 'TOP         = ', ma%top
      write(ndebug+myid,*) 'TOPRIGHT    = ', ma%topright
      write(ndebug+myid,*) 'RIGHT       = ', ma%right
      write(ndebug+myid,*) 'BOTTOMRIGHT = ', ma%bottomright
      write(ndebug+myid,*) 'BOTTOM      = ', ma%bottom
      write(ndebug+myid,*) 'BOTTOMLEFT  = ', ma%bottomleft
      write(ndebug+myid,*) 'LEFT        = ', ma%left
      write(ndebug+myid,*) 'GLOBAL J = ', global_dot_jstart , global_dot_jend
      write(ndebug+myid,*) 'GLOBAL I = ', global_dot_istart , global_dot_iend
      write(ndebug+myid,*) 'DOTPEXTJI1 : ', jde1 , jde2 , ide1 , ide2
      write(ndebug+myid,*) 'DOTPINTJI1 : ', jdi1 , jdi2 , idi1 , idi2
      write(ndebug+myid,*) 'DOTPINTJI2 : ', jdii1 , jdii2 , idii1 , idii2
      write(ndebug+myid,*) 'CRXPEXTJI1 : ', jce1 , jce2 , ice1 , ice2
      write(ndebug+myid,*) 'CRXPINTJI1 : ', jci1 , jci2 , ici1 , ici2
      write(ndebug+myid,*) 'CRXPINTJI2 : ', jcii1 , jcii2 , icii1 , icii2
      write(ndebug+myid,*) 'TOPBDY   : ', ma%has_bdytop
      write(ndebug+myid,*) 'BTMBDY   : ', ma%has_bdybottom
      write(ndebug+myid,*) 'RGTBDY   : ', ma%has_bdyright
      write(ndebug+myid,*) 'LFTBDY   : ', ma%has_bdyleft
      flush(ndebug+myid)
#endif
    end subroutine setup_model_indexes

    subroutine setup_boundaries(ldot,ba)
      implicit none
      logical , intent(in) :: ldot
      type(bound_area) , intent(out) :: ba
      integer(ik4) :: ic
      integer(ik4) :: igbb1 , igbb2 , igbt1 , igbt2
      integer(ik4) :: jgbl1 , jgbl2 , jgbr1 , jgbr2
      integer(ik4) :: i , j , i1 , i2 , j1 , j2 , iglob , jglob

      ba%dotflag = ldot
      call getmem2d(ba%ibnd,jde1,jde2,ide1,ide2,'setup_boundaries:ibnd')
      call getmem2d(ba%bsouth,jde1,jde2,ide1,ide2,'setup_boundaries:bsouth')
      call getmem2d(ba%bnorth,jde1,jde2,ide1,ide2,'setup_boundaries:bnorth')
      call getmem2d(ba%beast,jde1,jde2,ide1,ide2,'setup_boundaries:beast')
      call getmem2d(ba%bwest,jde1,jde2,ide1,ide2,'setup_boundaries:bwest')
      if ( ldot ) then
        ic = 0
        ba%nsp = nspgd
      else
        ic = 1
        ba%nsp = nspgx
      end if
      ba%ibnd(:,:) = -1
      igbb1 = 2
      igbb2 = ba%nsp-1
      jgbl1 = 2
      jgbl2 = ba%nsp-1
      igbt1 = iy-ic-ba%nsp+2
      igbt2 = iy-ic-1
      jgbr1 = jx-ic-ba%nsp+2
      jgbr2 = jx-ic-1
      i1 = 1
      i2 = iyp
      j1 = 1
      j2 = jxp
      ba%ns = 0
      ba%nn = 0
      ba%nw = 0
      ba%ne = 0
      ba%bsouth(:,:) = .false.
      ba%bnorth(:,:) = .false.
      ba%bwest(:,:) = .false.
      ba%beast(:,:) = .false.
      if ( ma%bandflag ) then
        ! Check for South boundary
        do i = i1 , i2
          iglob = global_dot_istart+i-1
          if ( iglob >= igbb1 .and. iglob <= igbb2 ) then
            do j = j1 , j2
              jglob = global_dot_jstart+j-1
              if ( jglob < jgbl1 .and. jglob > jgbr2 ) cycle
              ba%ibnd(j,i) = iglob-igbb1+2
              ba%bsouth(j,i) = .true.
              ba%ns = ba%ns+1
            end do
          end if
        end do
        ! North Boundary
        do i = i1 , i2
          iglob = global_dot_istart+i-1
          if ( iglob >= igbt1 .and. iglob <= igbt2 ) then
            do j = j1 , j2
              jglob = global_dot_jstart+j-1
              if ( jglob < jgbl1 .and. jglob > jgbr2 ) cycle
              ba%ibnd(j,i) = igbt2-iglob+2
              ba%bnorth(j,i) = .true.
              ba%nn = ba%nn+1
            end do
          end if
        end do
      else
        ! Check for South boundary
        do i = i1 , i2
          iglob = global_dot_istart+i-1
          if ( iglob >= igbb1 .and. iglob <= igbb2 ) then
            do j = j1 , j2
              jglob = global_dot_jstart+j-1
              if ( jglob >= jgbl1 .and. jglob <= jgbr2 ) then
                if ( jglob <= jgbl2 .and. iglob >= jglob ) cycle
                if ( jglob >= jgbr1 .and. iglob >= (jgbr2-jglob+2) ) cycle
                ba%ibnd(j,i) = iglob-igbb1+2
                ba%bsouth(j,i) = .true.
                ba%ns = ba%ns+1
              end if
            end do
          end if
        end do
        ! North Boundary
        do i = i1 , i2
          iglob = global_dot_istart+i-1
          if ( iglob >= igbt1 .and. iglob <= igbt2 ) then
            do j = j1 , j2
              jglob = global_dot_jstart+j-1
              if ( jglob >= jgbl1 .and. jglob <= jgbr2 ) then
                if ( jglob <= jgbl2 .and. (igbt2-iglob+2) >= jglob ) cycle
                if ( jglob >= jgbr1 .and. (igbt2-iglob) >= (jgbr2-jglob) ) cycle
                ba%ibnd(j,i) = igbt2-iglob+2
                ba%bnorth(j,i) = .true.
                ba%nn = ba%nn+1
              end if
            end do
          end if
        end do
        ! West boundary
        do i = i1 , i2
          iglob = global_dot_istart+i-1
          if ( iglob < igbb1 .or. iglob > igbt2 ) cycle
          do j = j1 , j2
            jglob = global_dot_jstart+j-1
            if ( jglob >= jgbl1 .and. jglob <= jgbl2 ) then
              if ( iglob < igbb2 .and. jglob > iglob ) cycle
              if ( iglob > igbt1 .and. jglob > (igbt2-iglob+2) ) cycle
              ba%ibnd(j,i) = jglob-jgbl1+2
              ba%bwest(j,i) = .true.
              ba%nw = ba%nw+1
            end if
          end do
        end do
        ! East boundary
        do i = i1 , i2
          iglob = global_dot_istart+i-1
          if ( iglob < igbb1 .or. iglob > igbt2 ) cycle
          do j = j1 , j2
            jglob = global_dot_jstart+j-1
            if ( jglob >= jgbr1 .and. jglob <= jgbr2 ) then
              if ( iglob < igbb2 .and. (jgbr2-jglob+2) > iglob ) cycle
              if ( iglob > igbt1 .and. (jgbr2-jglob) > (igbt2-iglob) ) cycle
              ba%ibnd(j,i) = jgbr2-jglob+2
              ba%beast(j,i) = .true.
              ba%ne = ba%ne+1
            end if
          end do
        end do
      end if
      ba%havebound = (ba%ns /= 0 .or. ba%nn /= 0 .or. &
                      ba%nw /= 0 .or. ba%ne /= 0)
#ifdef DEBUG
      write(ndebug+myid,*) 'DOT  : ', ldot
      write(ndebug+myid,*) 'BDYS : ', ba%ns
      write(ndebug+myid,*) 'BDYN : ', ba%nn
      write(ndebug+myid,*) 'BDYW : ', ba%nw
      write(ndebug+myid,*) 'BDYE : ', ba%ne

      do i = ide2 , ide1 , -1
        do j = jde1, jde2
          if ( ba%bsouth(j,i) ) then
            write(ndebug+myid,'(1a,i0.4)',advance='no') 'S' , ba%ibnd(j,i)
          else if ( ba%bnorth(j,i) ) then
            write(ndebug+myid,'(1a,i0.4)',advance='no') 'N' , ba%ibnd(j,i)
          else if ( ba%bwest(j,i) ) then
            write(ndebug+myid,'(1a,i0.4)',advance='no') 'W' , ba%ibnd(j,i)
          else if ( ba%beast(j,i) ) then
            write(ndebug+myid,'(1a,i0.4)',advance='no') 'E' , ba%ibnd(j,i)
          else
            write(ndebug+myid,'(1a,i0.4)',advance='no') 'X', 0
          end if
        end do
        write(ndebug+myid,*) ' '
      end do
#endif
    end subroutine setup_boundaries

    subroutine allocate_v3dbound(xb,ke,ldot)
      implicit none
      type(v3dbound) , intent(out) :: xb
      integer(ik4) , intent(in) :: ke
      logical , intent(in) :: ldot
      if ( ldot ) then
        call getmem3d(xb%b0,jde1-ma%jbl1,jde2+ma%jbr1, &
                            ide1-ma%ibb1,ide2+ma%ibt1,1,ke,'v3dbound:b0')
        call getmem3d(xb%b1,jde1-ma%jbl1,jde2+ma%jbr1, &
                            ide1-ma%ibb1,ide2+ma%ibt1,1,ke,'v3dbound:b1')
        call getmem3d(xb%bt,jde1-ma%jbl1,jde2+ma%jbr1, &
                            ide1-ma%ibb1,ide2+ma%ibt1,1,ke,'v3dbound:bt')
      else
        call getmem3d(xb%b0,jce1-ma%jbl1,jce2+ma%jbr1, &
                            ice1-ma%ibb1,ice2+ma%ibt1,1,ke,'v3dbound:b0')
        call getmem3d(xb%b1,jce1-ma%jbl1,jce2+ma%jbr1, &
                            ice1-ma%ibb1,ice2+ma%ibt1,1,ke,'v3dbound:b1')
        call getmem3d(xb%bt,jce1-ma%jbl1,jce2+ma%jbr1, &
                            ice1-ma%ibb1,ice2+ma%ibt1,1,ke,'v3dbound:bt')
      end if
    end subroutine allocate_v3dbound
!
    subroutine allocate_v2dbound(xb,ldot)
      implicit none
      type(v2dbound) , intent(out) :: xb
      logical , intent(in) :: ldot
      if ( ldot ) then
        call getmem2d(xb%b0,jde1-ma%jbl1,jde2+ma%jbr1, &
                            ide1-ma%ibb1,ide2+ma%ibt1,'v2dbound:b0')
        call getmem2d(xb%b1,jde1-ma%jbl1,jde2+ma%jbr1, &
                            ide1-ma%ibb1,ide2+ma%ibt1,'v2dbound:b1')
        call getmem2d(xb%bt,jde1-ma%jbl1,jde2+ma%jbr1, &
                            ide1-ma%ibb1,ide2+ma%ibt1,'v2dbound:bt')
      else
        call getmem2d(xb%b0,jce1-ma%jbl1,jce2+ma%jbr1, &
                            ice1-ma%ibb1,ice2+ma%ibt1,'v2dbound:b0')
        call getmem2d(xb%b1,jce1-ma%jbl1,jce2+ma%jbr1, &
                            ice1-ma%ibb1,ice2+ma%ibt1,'v2dbound:b1')
        call getmem2d(xb%bt,jce1-ma%jbl1,jce2+ma%jbr1, &
                            ice1-ma%ibb1,ice2+ma%ibt1,'v2dbound:bt')
      end if
    end subroutine allocate_v2dbound

    subroutine allocate_atmstate_a(atm,exchange_points)
      implicit none
      integer(ik4) , intent(in) :: exchange_points
      type(atmstate_a) , intent(out) :: atm
      integer(ik4) :: ib , it , jr , jl
      if ( exchange_points == zero_exchange_point ) then
        ib = 0
        it = 0
        jl = 0
        jr = 0
      else if ( exchange_points == one_exchange_point ) then
        ib = ma%ibb1
        it = ma%ibt1
        jl = ma%jbl1
        jr = ma%jbr1
      else if ( exchange_points == two_exchange_point ) then
        ib = ma%ibb2
        it = ma%ibt2
        jl = ma%jbl2
        jr = ma%jbr2
      else if ( exchange_points == four_exchange_point ) then
        ib = ma%ibb4
        it = ma%ibt4
        jl = ma%jbl4
        jr = ma%jbr4
      else
        ib = -1
        it = -1
        jl = -1
        jr = -1
        call fatal(__FILE__,__LINE__,'Uncoded number of exchange points')
      end if
      call getmem3d(atm%u,jde1-jl,jde2+jr,ide1-ib,ide2+it,1,kz,'atmstate:u')
      call getmem3d(atm%v,jde1-jl,jde2+jr,ide1-ib,ide2+it,1,kz,'atmstate:v')
      call getmem3d(atm%t,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:t')
      call getmem4d(atm%qx,jce1-jl,jce2+jr, &
                           ice1-ib,ice2+it,1,kz,1,nqx,'atmstate:qx')
      if ( ibltyp == 2 ) then
        call getmem3d(atm%tke,jce1-jl,jce2+jr,ice1-ib,ice2+it, &
                      1,kzp1,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem3d(atm%pr,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1-jl,jce2+jr, &
                ice1-ib,ice2+it,1,kz,'atmstate:rho')
        call getmem3d(atm%pp,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:pp')
        call getmem3d(atm%w,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kzp1,'atmstate:w')
      else
        call getmem3d(atm%pr,jce1,jce2,ice1,ice2,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1,jce2,ice1,ice2,1,kz,'atmstate:rho')
      end if
    end subroutine allocate_atmstate_a

    subroutine allocate_atmstate_b(atm,exchange_points)
      implicit none
      integer(ik4) , intent(in) :: exchange_points
      type(atmstate_b) , intent(out) :: atm
      integer(ik4) :: ib , it , jr , jl
      if ( exchange_points == zero_exchange_point ) then
        ib = 0
        it = 0
        jl = 0
        jr = 0
      else if ( exchange_points == one_exchange_point ) then
        ib = ma%ibb1
        it = ma%ibt1
        jl = ma%jbl1
        jr = ma%jbr1
      else if ( exchange_points == two_exchange_point ) then
        ib = ma%ibb2
        it = ma%ibt2
        jl = ma%jbl2
        jr = ma%jbr2
      else if ( exchange_points == four_exchange_point ) then
        ib = ma%ibb4
        it = ma%ibt4
        jl = ma%jbl4
        jr = ma%jbr4
      else
        ib = -1
        it = -1
        jl = -1
        jr = -1
        call fatal(__FILE__,__LINE__,'Uncoded number of exchange points')
      end if
      call getmem3d(atm%u,jde1-jl,jde2+jr,ide1-ib,ide2+it,1,kz,'atmstate:u')
      call getmem3d(atm%v,jde1-jl,jde2+jr,ide1-ib,ide2+it,1,kz,'atmstate:v')
      call getmem3d(atm%t,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:t')
      call getmem4d(atm%qx,jce1-jl,jce2+jr, &
                           ice1-ib,ice2+it,1,kz,1,nqx,'atmstate:qx')
      if ( ibltyp == 2 ) then
        call getmem3d(atm%tke,jce1-jl,jce2+jr,ice1-ib,ice2+it, &
                      1,kzp1,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem3d(atm%pr,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1-jl,jce2+jr, &
                ice1-ib,ice2+it,1,kz,'atmstate:rho')
        call getmem3d(atm%pp,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:pp')
        call getmem3d(atm%w,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kzp1,'atmstate:w')
      else
        call getmem3d(atm%pr,jce1,jce2,ice1,ice2,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1,jce2,ice1,ice2,1,kz,'atmstate:rho')
      end if
    end subroutine allocate_atmstate_b

    subroutine allocate_atmstate_c(atm,exchange_points)
      implicit none
      integer(ik4) , intent(in) :: exchange_points
      type(atmstate_c) , intent(out) :: atm
      integer(ik4) :: ib , it , jr , jl
      if ( exchange_points == zero_exchange_point ) then
        ib = 0
        it = 0
        jl = 0
        jr = 0
      else if ( exchange_points == one_exchange_point ) then
        ib = ma%ibb1
        it = ma%ibt1
        jl = ma%jbl1
        jr = ma%jbr1
      else if ( exchange_points == two_exchange_point ) then
        ib = ma%ibb2
        it = ma%ibt2
        jl = ma%jbl2
        jr = ma%jbr2
      else if ( exchange_points == four_exchange_point ) then
        ib = ma%ibb4
        it = ma%ibt4
        jl = ma%jbl4
        jr = ma%jbr4
      else
        ib = -1
        it = -1
        jl = -1
        jr = -1
        call fatal(__FILE__,__LINE__,'Uncoded number of exchange points')
      end if
      call getmem3d(atm%u,jde1-jl,jde2+jr,ide1-ib,ide2+it,1,kz,'atmstate:u')
      call getmem3d(atm%v,jde1-jl,jde2+jr,ide1-ib,ide2+it,1,kz,'atmstate:v')
      call getmem3d(atm%t,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:t')
      call getmem4d(atm%qx,jce1-jl,jce2+jr, &
                           ice1-ib,ice2+it,1,kz,1,nqx,'atmstate:qx')
      if ( ibltyp == 2 ) then
        call getmem3d(atm%tke,jce1-jl,jce2+jr,ice1-ib,ice2+it, &
                      1,kzp1,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem3d(atm%pr,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1-jl,jce2+jr, &
                ice1-ib,ice2+it,1,kz,'atmstate:rho')
        call getmem3d(atm%pp,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:pp')
        call getmem3d(atm%w,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kzp1,'atmstate:w')
      else
        call getmem3d(atm%pr,jce1,jce2,ice1,ice2,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1,jce2,ice1,ice2,1,kz,'atmstate:rho')
      end if
    end subroutine allocate_atmstate_c

    subroutine allocate_atmstate_decoupled(atm,exchange_points)
      implicit none
      integer(ik4) , intent(in) :: exchange_points
      type(atmstate_decoupled) , intent(out) :: atm
      integer(ik4) :: ib , it , jr , jl
      if ( exchange_points == zero_exchange_point ) then
        ib = 0
        it = 0
        jl = 0
        jr = 0
      else if ( exchange_points == one_exchange_point ) then
        ib = ma%ibb1
        it = ma%ibt1
        jl = ma%jbl1
        jr = ma%jbr1
      else if ( exchange_points == two_exchange_point ) then
        ib = ma%ibb2
        it = ma%ibt2
        jl = ma%jbl2
        jr = ma%jbr2
      else if ( exchange_points == four_exchange_point ) then
        ib = ma%ibb4
        it = ma%ibt4
        jl = ma%jbl4
        jr = ma%jbr4
      else
        ib = -1
        it = -1
        jl = -1
        jr = -1
        call fatal(__FILE__,__LINE__,'Uncoded number of exchange points')
      end if
      call getmem3d(atm%u,jde1-jl,jde2+jr,ide1-ib,ide2+it,1,kz,'atmstate:u')
      call getmem3d(atm%v,jde1-jl,jde2+jr,ide1-ib,ide2+it,1,kz,'atmstate:v')
      call getmem3d(atm%t,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:t')
      call getmem4d(atm%qx,jce1-jl,jce2+jr, &
                           ice1-ib,ice2+it,1,kz,1,nqx,'atmstate:qx')
      if ( ibltyp == 2 ) then
        call getmem3d(atm%tke,jce1-jl,jce2+jr,ice1-ib,ice2+it, &
                      1,kzp1,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem3d(atm%pr,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1-jl,jce2+jr, &
                ice1-ib,ice2+it,1,kz,'atmstate:rho')
        call getmem3d(atm%pp,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:pp')
        call getmem3d(atm%w,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kzp1,'atmstate:w')
      else
        call getmem3d(atm%pr,jce1,jce2,ice1,ice2,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1,jce2,ice1,ice2,1,kz,'atmstate:rho')
      end if
    end subroutine allocate_atmstate_decoupled

    subroutine allocate_atmstate_tendency(atm,exchange_points)
      implicit none
      integer(ik4) , intent(in) :: exchange_points
      type(atmstate_tendency) , intent(out) :: atm
      integer(ik4) :: ib , it , jr , jl
      if ( exchange_points == zero_exchange_point ) then
        ib = 0
        it = 0
        jl = 0
        jr = 0
      else if ( exchange_points == one_exchange_point ) then
        ib = ma%ibb1
        it = ma%ibt1
        jl = ma%jbl1
        jr = ma%jbr1
      else if ( exchange_points == two_exchange_point ) then
        ib = ma%ibb2
        it = ma%ibt2
        jl = ma%jbl2
        jr = ma%jbr2
      else if ( exchange_points == four_exchange_point ) then
        ib = ma%ibb4
        it = ma%ibt4
        jl = ma%jbl4
        jr = ma%jbr4
      else
        ib = -1
        it = -1
        jl = -1
        jr = -1
        call fatal(__FILE__,__LINE__,'Uncoded number of exchange points')
      end if
      call getmem3d(atm%u,jde1-jl,jde2+jr,ide1-ib,ide2+it,1,kz,'atmstate:u')
      call getmem3d(atm%v,jde1-jl,jde2+jr,ide1-ib,ide2+it,1,kz,'atmstate:v')
      call getmem3d(atm%t,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:t')
      call getmem4d(atm%qx,jce1-jl,jce2+jr, &
                           ice1-ib,ice2+it,1,kz,1,nqx,'atmstate:qx')
      if ( ibltyp == 2 ) then
        call getmem3d(atm%tke,jce1-jl,jce2+jr,ice1-ib,ice2+it, &
                      1,kzp1,'atmstate:tke')
      end if
      if ( idynamic == 2 ) then
        call getmem3d(atm%pr,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1-jl,jce2+jr, &
                ice1-ib,ice2+it,1,kz,'atmstate:rho')
        call getmem3d(atm%pp,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:pp')
        call getmem3d(atm%w,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kzp1,'atmstate:w')
      else
        call getmem3d(atm%pr,jce1,jce2,ice1,ice2,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1,jce2,ice1,ice2,1,kz,'atmstate:rho')
      end if
    end subroutine allocate_atmstate_tendency

    subroutine allocate_reference_atmosphere(atm)
      implicit none
      type(reference_atmosphere) , intent(out) :: atm
        call getmem3d(atm%t,jce1,jce2,ice1,ice2,1,kz,'atmstate:t')
        call getmem3d(atm%pr,jce1,jce2,ice1,ice2,1,kz,'atmstate:pr')
        call getmem3d(atm%rho,jce1,jce2,ice1,ice2,1,kz,'atmstate:rho')
    end subroutine allocate_reference_atmosphere

    subroutine allocate_tendiag(dia)
      implicit none
      type(tendiag) , intent(out) :: dia
      call getmem3d(dia%adh,jce1,jce2,ice1,ice2,1,kz,'tendiag:adh')
      call getmem3d(dia%adv,jce1,jce2,ice1,ice2,1,kz,'tendiag:adv')
      call getmem3d(dia%tbl,jce1,jce2,ice1,ice2,1,kz,'tendiag:tbl')
      call getmem3d(dia%con,jce1,jce2,ice1,ice2,1,kz,'tendiag:con')
      call getmem3d(dia%bdy,jce1,jce2,ice1,ice2,1,kz,'tendiag:bdy')
      call getmem3d(dia%adi,jce1,jce2,ice1,ice2,1,kz,'tendiag:adi')
      call getmem3d(dia%dif,jce1,jce2,ice1,ice2,1,kz,'tendiag:dif')
      call getmem3d(dia%rad,jce1,jce2,ice1,ice2,1,kz,'tendiag:rad')
      call getmem3d(dia%lsc,jce1,jce2,ice1,ice2,1,kz,'tendiag:lsc')
    end subroutine allocate_tendiag

    subroutine allocate_qendiag(dia)
      implicit none
      type(qendiag) , intent(out) :: dia
      call getmem3d(dia%adh,jce1,jce2,ice1,ice2,1,kz,'tendiag:adh')
      call getmem3d(dia%adv,jce1,jce2,ice1,ice2,1,kz,'tendiag:adv')
      call getmem3d(dia%tbl,jce1,jce2,ice1,ice2,1,kz,'tendiag:tbl')
      call getmem3d(dia%con,jce1,jce2,ice1,ice2,1,kz,'tendiag:con')
      call getmem3d(dia%bdy,jce1,jce2,ice1,ice2,1,kz,'tendiag:bdy')
      call getmem3d(dia%adi,jce1,jce2,ice1,ice2,1,kz,'tendiag:adi')
      call getmem3d(dia%dif,jce1,jce2,ice1,ice2,1,kz,'tendiag:dif')
      call getmem3d(dia%rad,jce1,jce2,ice1,ice2,1,kz,'tendiag:rad')
      call getmem3d(dia%lsc,jce1,jce2,ice1,ice2,1,kz,'tendiag:lsc')
    end subroutine allocate_qendiag

    subroutine allocate_domain(dom)
      implicit none
      type(domain) , intent(out) :: dom
      call getmem2d(dom%ht,jde1-ma%jbl1,jde2+ma%jbr1, &
                           ide1-ma%ibb1,ide2+ma%ibt1,'storage:ht')
      call getmem2d(dom%lndcat,jde1,jde2,ide1,ide2,'storage:lndcat')
      call getmem2d(dom%xlat,jde1,jde2,ide1,ide2,'storage:xlat')
      call getmem2d(dom%xlon,jde1,jde2,ide1,ide2,'storage:xlon')
      call getmem2d(dom%mask,jde1,jde2,ide1,ide2,'storage:mask')
      call getmem2d(dom%dlat,jde1,jde2,ide1,ide2,'storage:dlat')
      call getmem2d(dom%dlon,jde1,jde2,ide1,ide2,'storage:dlon')
      call getmem2d(dom%msfx,jde1-ma%jbl2,jde2+ma%jbr2, &
                             ide1-ma%ibb2,ide2+ma%ibt2,'storage:msfx')
      call getmem2d(dom%msfd,jde1-ma%jbl2,jde2+ma%jbr2, &
                             ide1-ma%ibb2,ide2+ma%ibt2,'storage:msfd')
      call getmem2d(dom%coriol,jde1,jde2,ide1,ide2,'storage:f')
      call getmem2d(dom%snowam,jde1,jde2,ide1,ide2,'storage:snowam')
      call getmem2d(dom%ldmsk,jci1,jci2,ici1,ici2,'storage:ldmsk')
      call getmem2d(dom%iveg,jci1,jci2,ici1,ici2,'storage:iveg')
      if ( lakemod == 1 ) then
        call getmem2d(dom%dhlake,jde1,jde2,ide1,ide2,'storage:dhlake')
      end if
    end subroutine allocate_domain
!
    subroutine allocate_domain_subgrid(sub)
      implicit none
      type(domain_subgrid) , intent(out) :: sub
      call getmem3d(sub%ht,1,nnsg,jde1,jde2,ide1,ide2,'storage:ht')
      call getmem3d(sub%lndcat,1,nnsg,jde1,jde2,ide1,ide2,'storage:lndcat')
      call getmem3d(sub%xlat,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlat')
      call getmem3d(sub%xlon,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlon')
      call getmem3d(sub%mask,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlon')
      call getmem3d(sub%ldmsk,1,nnsg,jci1,jci2,ici1,ici2,'storage:ldmsk')
      call getmem3d(sub%iveg,1,nnsg,jci1,jci2,ici1,ici2,'storage:iveg')
      if ( lakemod == 1 ) then
        call getmem3d(sub%dhlake,1,nnsg,jde1,jde2,ide1,ide2,'storage:dhlake')
      end if
    end subroutine allocate_domain_subgrid

    subroutine allocate_surfstate(sfs)
      implicit none
      type(surfstate) , intent(out) :: sfs
      call getmem2d(sfs%psa,jce1-ma%jbl1,jce2+ma%jbr1, &
                            ice1-ma%ibb1,ice2+ma%ibt1,'surf:psa')
      call getmem2d(sfs%psb,jce1-ma%jbl1,jce2+ma%jbr1, &
                            ice1-ma%ibb1,ice2+ma%ibt1,'surf:psb')
      call getmem2d(sfs%psc,jce1,jce2,ice1,ice2,'surf:psc')
      call getmem2d(sfs%psdota,jde1-ma%jbl1,jde2+ma%jbr1, &
                            ide1-ma%ibb1,ide2+ma%ibt1,'surf:psdota')
      call getmem2d(sfs%psdotb,jde1-ma%jbl1,jde2+ma%jbr1, &
                            ide1-ma%ibb1,ide2+ma%ibt1,'surf:psdotb')
      call getmem2d(sfs%tga,jci1,jci2,ici1,ici2,'surf:tga')
      call getmem2d(sfs%tgb,jci1,jci2,ici1,ici2,'surf:tgb')
      call getmem2d(sfs%hfx,jci1,jci2,ici1,ici2,'surf:hfx')
      call getmem2d(sfs%qfx,jci1,jci2,ici1,ici2,'surf:qfx')
      call getmem2d(sfs%rainc,jci1,jci2,ici1,ici2,'surf:rainc')
      call getmem2d(sfs%rainnc,jci1,jci2,ici1,ici2,'surf:rainnc')
      if ( ipptls == 2 ) then
        call getmem2d(sfs%snownc,jci1,jci2,ici1,ici2,'surf:snownc')
      end if
      call getmem2d(sfs%tgbb,jci1,jci2,ici1,ici2,'surf:tgbb')
      call getmem2d(sfs%uvdrag,jci1,jci2,ici1,ici2,'surf:uvdrag')
    end subroutine allocate_surfstate
!
    subroutine allocate_slice(ax)
      implicit none
      type(slice) , intent(out) :: ax
      call getmem3d(ax%pf3d,jce1,jce2,ice1,ice2,1,kz+1,'slice:pb3d')
      call getmem3d(ax%qsb3d,jce1,jce2,ice1,ice2,1,kz,'slice:qsb3d')
      call getmem3d(ax%rhb3d,jce1,jce2,ice1,ice2,1,kz,'slice:rhb3d')
      call getmem3d(ax%thx3d,jce1,jce2,ice1,ice2,1,kz,'slice:thx3d')
      call getmem3d(ax%ubx3d,jce1-ma%jbl2,jce2+ma%jbr2, &
                             ice1-ma%ibb2,ice2+ma%ibt2,1,kz,'slice:ubx3d')
      call getmem3d(ax%vbx3d,jce1-ma%jbl2,jce2+ma%jbr2, &
                             ice1-ma%ibb2,ice2+ma%ibt2,1,kz,'slice:vbx3d')
      call getmem4d(ax%qxb3d,jce1-ma%jbl2,jce2+ma%jbr2, &
                             ice1-ma%ibb2,ice2+ma%ibt2,1,kz,1,nqx,'slice:qxb3d')
      call getmem3d(ax%tb3d,jce1-ma%jbl2,jce2+ma%jbr2, &
                            ice1-ma%ibb2,ice2+ma%ibt2,1,kz,'slice:tb3d')
      call getmem3d(ax%ubd3d,jde1-ma%jbl2,jde2+ma%jbr2, &
                             ide1-ma%ibb2,ide2+ma%ibt2,1,kz,'slice:ubd3d')
      call getmem3d(ax%vbd3d,jde1-ma%jbl2,jde2+ma%jbr2, &
                             ide1-ma%ibb2,ide2+ma%ibt2,1,kz,'slice:vbd3d')
      call getmem3d(ax%zq,jce1,jce2,ice1,ice2,1,kzp1,'slice:zq')
      call getmem3d(ax%za,jce1,jce2,ice1,ice2,1,kz,'slice:za')
      call getmem3d(ax%dzq,jce1,jce2,ice1,ice2,1,kz,'slice:dzq')
      call getmem2d(ax%rhox2d,jce1,jce2,ice1,ice2,'slice:rhox2d')
      call getmem2d(ax%ps2d,jce1,jce2,ice1,ice2,'slice:ps2d')
      call getmem3d(ax%wpx3d,jci1,jci2,ici1,ici2,1,kz,'slice:wpx3d')
      if ( ichem == 1 ) then
        call getmem4d(ax%chib3d,jce1-ma%jbl2,jce2+ma%jbr2, &
                                ice1-ma%ibb2,ice2+ma%ibt2, &
                                1,kz,1,ntr,'slice:chib3d')
      end if
      if ( ibltyp == 2 ) then
        call getmem3d(ax%tkeb3d,jce1-ma%jbl2,jce2+ma%jbr2, &
                                ice1-ma%ibb2,ice2+ma%ibt2,1,kzp1,'slice:tkeb3d')
      end if
    end subroutine allocate_slice

    subroutine allocate_diffx(dx)
      implicit none
      type(diffx) , intent(out) :: dx
      call getmem3d(dx%difuu,jdi1,jdi2,idi1,idi2,1,kz,'diffx:difuu')
      call getmem3d(dx%difuv,jdi1,jdi2,idi1,idi2,1,kz,'diffx:difuv')
      call getmem3d(dx%difft,jci1,jci2,ici1,ici2,1,kz,'diffx:difft')
      call getmem4d(dx%diffqx,jci1,jci2,ici1,ici2,1,kz,1,nqx,'diffx:diffqx')
    end subroutine allocate_diffx

    subroutine allocate_mod_atm_interface
      implicit none

      call allocate_domain(mddom)
      call allocate_domain_subgrid(mdsub)

      if ( isladvec == 1 ) then
        call allocate_atmstate_a(atm1,four_exchange_point)
        call allocate_atmstate_b(atm2,four_exchange_point)
        call allocate_atmstate_decoupled(atmx,four_exchange_point)
      else
        call allocate_atmstate_a(atm1,one_exchange_point)
        call allocate_atmstate_b(atm2,one_exchange_point)
        call allocate_atmstate_decoupled(atmx,one_exchange_point)
      end if
      call allocate_atmstate_c(atmc,zero_exchange_point)
      call allocate_atmstate_tendency(aten,zero_exchange_point)
      if ( ibltyp == 2 ) then
        call allocate_atmstate_tendency(uwten,one_exchange_point)
      end if

      call allocate_surfstate(sfs)

      call allocate_slice(atms)

      call allocate_diffx(adf)

      ! FAB:
      !    complete for diag on water quantitiies idiag = 2, 3 etc
      if ( idiag > 0 ) then
        call allocate_tendiag(tdiag)
        call allocate_qendiag(qdiag)
      end if

      call getmem2d(ts0,jce1,jce2,ice1,ice2,'storage:ts0')
      call getmem2d(ts1,jce1,jce2,ice1,ice2,'storage:ts1')

      call getmem3d(dstor,jde1,jde2,ide1,ide2,1,nsplit,'storage:dstor')
      call getmem3d(hstor,jde1,jde2,ide1,ide2,1,nsplit,'storage:hstor')

      call getmem2d(hgfact,jde1,jde2,ide1,ide2,'storage:hgfact')
      call getmem3d(omega,jci1,jci2,ici1,ici2,1,kz,'storage:omega')
      call getmem3d(qdot,jce1-ma%jbl1,jce2+ma%jbr1, &
                         ice1-ma%ibb1,ice2+ma%ibt1,1,kzp1,'storage:qdot')
      call getmem2d(ktrop,jci1,jci2,ici1,ici2,'storage:ktrop')
      call getmem2d(coszrs,jci1,jci2,ici1,ici2,'storage:coszrs')
      call getmem2d(pptc,jci1,jci2,ici1,ici2,'storage:pptc')
      call getmem2d(prca,jci1,jci2,ici1,ici2,'storage:prca')
      call getmem2d(icumbot,jci1,jci2,ici1,ici2,'storage:icumbot')
      call getmem2d(icumtop,jci1,jci2,ici1,ici2,'storage:icumtop')
      if ( ichem == 1 ) then
        call getmem3d(convpr,jci1,jci2,ici1,ici2,1,kz,'storage:convpr')
        call getmem2d(ssw2da,jci1,jci2,ici1,ici2,'storage:ssw2da')
        call getmem2d(sfracv2d,jci1,jci2,ici1,ici2,'storage:sfracv2d')
        call getmem2d(sfracb2d,jci1,jci2,ici1,ici2,'storage:sfracb2d')
        call getmem2d(sfracs2d,jci1,jci2,ici1,ici2,'storage:sfracs2d')
        call getmem2d(svegfrac2d,jci1,jci2,ici1,ici2,'storage:svegfrac2d')
        call getmem2d(sxlai2d,jci1,jci2,ici1,ici2,'storage:sxlai2d')
      end if
      if ( iocncpl == 1 ) then
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
      call getmem2d(flw,jci1,jci2,ici1,ici2,'storage:flw')
      call getmem2d(flwd,jci1,jci2,ici1,ici2,'storage:flwd')
      call getmem2d(fsw,jci1,jci2,ici1,ici2,'storage:fsw')
      call getmem2d(sabveg,jci1,jci2,ici1,ici2,'storage:sabveg')
      call getmem2d(solis,jci1,jci2,ici1,ici2,'storage:solis')
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

      call assignpnt(atm2%pr,atms%pb3d)
      call assignpnt(atm2%rho,atms%rhob3d)
      call getmem2d(zpbl,jci1,jci2,ici1,ici2,'storage:zpbl')
      call getmem2d(kpbl,jci1,jci2,ici1,ici2,'storage:kpbl')

      if ( idynamic == 2 ) then
        call allocate_reference_atmosphere(atm0)
        call getmem2d(dpsdxm,jci1,jci2,ici1,ici2,'storage:dpsdxm')
        call getmem2d(dpsdym,jci1,jci2,ici1,ici2,'storage:dpsdym')
      end if

    end subroutine allocate_mod_atm_interface

end module mod_atm_interface

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
