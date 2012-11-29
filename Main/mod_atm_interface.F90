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
!
  use mod_dynparam
  use mod_stdio
  use mod_constants , only : d_rfour
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_service
  use mod_memutil

  private

  logical , public , parameter :: cross = .false.
  logical , public , parameter :: dot = .true.
!
! Storage for all the 3d prognostic variables in two
! timesteps and all the 2d variables and constants
!
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

  type surfstate
    real(rk8) , pointer , dimension(:,:) :: psa
    real(rk8) , pointer , dimension(:,:) :: psb
    real(rk8) , pointer , dimension(:,:) :: tga
    real(rk8) , pointer , dimension(:,:) :: tgb
    real(rk8) , pointer , dimension(:,:) :: rainc
    real(rk8) , pointer , dimension(:,:) :: rainnc
    real(rk8) , pointer , dimension(:,:) :: hfx
    real(rk8) , pointer , dimension(:,:) :: qfx
    real(rk8) , pointer , dimension(:,:) :: tgbb
    real(rk8) , pointer , dimension(:,:) :: uvdrag
  end type surfstate

  type slice
    real(rk8) , pointer , dimension(:,:,:) :: tb3d
    real(rk8) , pointer , dimension(:,:,:) :: thx3d
    real(rk8) , pointer , dimension(:,:,:) :: pb3d
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

  public :: atmstate , domain , surfstate , slice
  public :: diffx , v2dbound , v3dbound , bound_area , model_area

  type(domain) , public :: mddom
  type(atmstate) , public :: atm1 , atm2
  type(atmstate) , public :: atmx , atmc , aten , holtten , uwten
  type(surfstate) , public :: sfs
  type(slice) , public :: atms
  type(diffx) , public :: adf
  type(v3dbound) , public :: xtb , xqb , xub , xvb
  type(v2dbound) , public :: xpsb
  type(bound_area) , public :: ba_cr , ba_dt

  public :: allocate_mod_atm_interface , allocate_atmstate , allocate_domain
  public :: allocate_surfstate , allocate_v3dbound , allocate_v2dbound
  public :: setup_boundaries , setup_model_indexes

  real(rk8) , public , pointer , dimension(:,:) :: hgfact
  real(rk8) , public , pointer , dimension(:,:) :: psdot
  real(rk8) , public , pointer , dimension(:,:,:) :: dstor
  real(rk8) , public , pointer , dimension(:,:,:) :: hstor
  real(rk8) , public , pointer , dimension(:,:) :: ts0 , ts1
!
  real(rk8) , public , pointer , dimension(:,:,:) :: qdot , omega
!
  integer(ik4) , public , parameter :: zero_exchange_point = 0
  integer(ik4) , public , parameter :: one_exchange_point = 1
  integer(ik4) , public , parameter :: two_exchange_point = 2

#ifdef DEBUG
  type(grid_nc_var4d) , public :: nc_4d
  type(grid_nc_var3d) , public :: nc_3d
  type(grid_nc_var2d) , public :: nc_2d
#endif

  contains 
!
    subroutine setup_model_indexes
      implicit none
      ma%jbl1 = 1
      ma%jbl2 = 2
      ma%jbr1 = 1
      ma%jbr2 = 2
      ma%ibt1 = 1
      ma%ibt2 = 2
      ma%ibb1 = 1
      ma%ibb2 = 2
      if ( ma%has_bdyleft ) then
        ma%jbl1 = 0
        ma%jbl2 = 0
      end if
      if ( ma%has_bdyright ) then
        ma%jbr1 = 0
        ma%jbr2 = 0
      end if
      if ( ma%has_bdytop ) then
        ma%ibt1 = 0
        ma%ibt2 = 0
      end if
      if ( ma%has_bdybottom ) then
        ma%ibb1 = 0
        ma%ibb2 = 0
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
!
    subroutine allocate_atmstate(atm,ibltyp,lpar,exchange_points)
      implicit none
      logical , intent(in) :: lpar
      integer(ik4) , intent(in) :: ibltyp
      integer(ik4) , intent(in) :: exchange_points
      type(atmstate) , intent(out) :: atm
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
      else
        call fatal(__FILE__,__LINE__,'Uncoded number of exchange points')
      end if

      if (lpar) then
        call getmem3d(atm%u,jde1-jl,jde2+jr,ide1-ib,ide2+it,1,kz,'atmstate:u')
        call getmem3d(atm%v,jde1-jl,jde2+jr,ide1-ib,ide2+it,1,kz,'atmstate:v')
        call getmem3d(atm%t,jce1-jl,jce2+jr,ice1-ib,ice2+it,1,kz,'atmstate:t')
        call getmem4d(atm%qx,jce1-jl,jce2+jr, &
                             ice1-ib,ice2+it,1,kz,1,nqx,'atmstate:qx')
        if ( ibltyp == 2 .or. ibltyp == 99 ) then
          call getmem3d(atm%tke,jce1-jl,jce2+jr,ice1-ib,ice2+it, &
                        1,kzp1,'atmstate:tke')
        end if
      else
        call getmem3d(atm%u,jdot1,jdot2,idot1,idot2,1,kz,'atmstate:u')
        call getmem3d(atm%v,jdot1,jdot2,idot1,idot2,1,kz,'atmstate:v')
        call getmem3d(atm%t,jcross1,jcross2,icross1,icross2,1,kz,'atmstate:t')
        call getmem4d(atm%qx,jcross1,jcross2, &
                             icross1,icross2,1,kz,1,nqx,'atmstate:qx')
        if ( ibltyp == 2 .or. ibltyp == 99 ) then
          call getmem3d(atm%tke, &
                        jcross1,jcross2,icross1,icross2,1,kzp1,'atmstate:tke')
        end if
      end if
    end subroutine allocate_atmstate
!
    subroutine allocate_domain(dom,lpar)
      implicit none
      logical , intent(in) :: lpar
      type(domain) , intent(out) :: dom

      if (lpar) then
        call getmem2d(dom%ht,jde1-ma%jbl1,jde2+ma%jbr1, &
                             ide1-ma%ibb1,ide2+ma%ibt1,'atm_interface:ht')
        call getmem2d(dom%lndcat,jde1,jde2,ide1,ide2,'atm_interface:lndcat')
        call getmem2d(dom%xlat,jde1,jde2,ide1,ide2,'atm_interface:xlat')
        call getmem2d(dom%xlon,jde1,jde2,ide1,ide2,'atm_interface:xlon')
        call getmem2d(dom%mask,jde1,jde2,ide1,ide2,'atm_interface:mask')
        call getmem2d(dom%dlat,jde1,jde2,ide1,ide2,'atm_interface:dlat')
        call getmem2d(dom%dlon,jde1,jde2,ide1,ide2,'atm_interface:dlon')
        call getmem2d(dom%msfx,jde1-ma%jbl2,jde2+ma%jbr2, &
                               ide1-ma%ibb2,ide2+ma%ibt2,'atm_interface:msfx')
        call getmem2d(dom%msfd,jde1-ma%jbl2,jde2+ma%jbr2, &
                               ide1-ma%ibb2,ide2+ma%ibt2,'atm_interface:msfd')
        call getmem2d(dom%coriol,jde1,jde2,ide1,ide2,'atm_interface:f')
        call getmem2d(dom%snowam,jde1,jde2,ide1,ide2,'atm_interface:snowam')
        if ( lakemod == 1 ) &
          call getmem2d(dom%dhlake,jde1,jde2,ide1,ide2,'atm_interface:dhlake')
      else
        call getmem2d(dom%ht,jdot1,jdot2,idot1,idot2,'atm_interface:ht')
        call getmem2d(dom%lndcat,jdot1,jdot2,idot1,idot2,'atm_interface:lndcat')
        call getmem2d(dom%xlat,jdot1,jdot2,idot1,idot2,'atm_interface:xlat')
        call getmem2d(dom%xlon,jdot1,jdot2,idot1,idot2,'atm_interface:xlon')
        call getmem2d(dom%mask,jdot1,jdot2,idot1,idot2,'atm_interface:mask')
        call getmem2d(dom%dlat,jdot1,jdot2,idot1,idot2,'atm_interface:dlat')
        call getmem2d(dom%dlon,jdot1,jdot2,idot1,idot2,'atm_interface:dlon')
        call getmem2d(dom%msfx,jdot1,jdot2,idot1,idot2,'atm_interface:msfx')
        call getmem2d(dom%msfd,jdot1,jdot2,idot1,idot2,'atm_interface:msfd')
        call getmem2d(dom%coriol,jdot1,jdot2,idot1,idot2,'atm_interface:f')
        call getmem2d(dom%snowam,jdot1,jdot2,idot1,idot2,'atm_interface:snowam')
        if ( lakemod == 1 ) then
          call getmem2d(dom%dhlake,jdot1,jdot2, &
            idot1,idot2,'atm_interface:dhlake')
        end if
      end if
    end subroutine allocate_domain
!
    subroutine allocate_surfstate(sfs,lpar)
      implicit none
      type(surfstate) , intent(out) :: sfs
      logical , intent(in) :: lpar
      if (lpar) then
        call getmem2d(sfs%psa,jce1-ma%jbl1,jce2+ma%jbr1, &
                              ice1-ma%ibb1,ice2+ma%ibt1,'surf:psa')
        call getmem2d(sfs%psb,jce1-ma%jbl1,jce2+ma%jbr1, &
                              ice1-ma%ibb1,ice2+ma%ibt1,'surf:psb')
        call getmem2d(sfs%tga,jce1,jce2,ice1,ice2,'surf:tga')
        call getmem2d(sfs%tgb,jce1,jce2,ice1,ice2,'surf:tgb')
        call getmem2d(sfs%hfx,jci1,jci2,ici1,ici2,'surf:hfx')
        call getmem2d(sfs%qfx,jci1,jci2,ici1,ici2,'surf:qfx')
        call getmem2d(sfs%rainc,jci1,jci2,ici1,ici2,'surf:rainc')
        call getmem2d(sfs%rainnc,jci1,jci2,ici1,ici2,'surf:rainnc')
        call getmem2d(sfs%tgbb,jci1,jci2,ici1,ici2,'surf:tgbb')
        call getmem2d(sfs%uvdrag,jci1,jci2,ici1,ici2,'surf:uvdrag')
      else
        call getmem2d(sfs%psa,jcross1,jcross2,icross1,icross2,'surf:psa')
        call getmem2d(sfs%psb,jcross1,jcross2,icross1,icross2,'surf:psb')
        call getmem2d(sfs%tga,jcross1,jcross2,icross1,icross2,'surf:tga')
        call getmem2d(sfs%tgb,jcross1,jcross2,icross1,icross2,'surf:tgb')
        call getmem2d(sfs%hfx,jcross1,jcross2,icross1,icross2,'surf:hfx')
        call getmem2d(sfs%qfx,jcross1,jcross2,icross1,icross2,'surf:qfx')
        call getmem2d(sfs%rainc,jcross1,jcross2,icross1,icross2,'surf:rainc')
        call getmem2d(sfs%rainnc,jcross1,jcross2,icross1,icross2,'surf:rainnc')
        call getmem2d(sfs%tgbb,jcross1,jcross2,icross1,icross2,'surf:tgbb')
        call getmem2d(sfs%uvdrag,jcross1,jcross2,icross1,icross2,'surf:uvdrag')
      end if
    end subroutine allocate_surfstate
!
    subroutine allocate_slice(ax,ibltyp)
      implicit none
      type(slice) , intent(out) :: ax
      integer(ik4) , intent(in) :: ibltyp
      call getmem3d(ax%pb3d,jce1,jce2,ice1,ice2,1,kz,'slice:pb3d')
      call getmem3d(ax%qsb3d,jce1,jce2,ice1,ice2,1,kz,'slice:qsb3d')
      call getmem3d(ax%rhb3d,jce1,jce2,ice1,ice2,1,kz,'slice:rhb3d')
      call getmem3d(ax%rhob3d,jce1,jce2,ice1,ice2,1,kz,'slice:rhob3d')
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
      if ( ichem == 1 ) then
        call getmem4d(ax%chib3d,jce1-ma%jbl2,jce2+ma%jbr2, &
                                ice1-ma%ibb2,ice2+ma%ibt2, &
                                1,kz,1,ntr,'slice:chib3d')
      end if
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call getmem3d(ax%tkeb3d,jce1-ma%jbl2,jce2+ma%jbr2, &
                                ice1-ma%ibb2,ice2+ma%ibt2,1,kzp1,'slice:tkeb3d')
      end if
    end subroutine allocate_slice
!
    subroutine allocate_diffx(dx)
      implicit none
      type(diffx) , intent(out) :: dx
      call getmem3d(dx%difuu,jdi1,jdi2,idi1,idi2,1,kz,'diffx:difuu')
      call getmem3d(dx%difuv,jdi1,jdi2,idi1,idi2,1,kz,'diffx:difuv')
      call getmem3d(dx%difft,jci1,jci2,ici1,ici2,1,kz,'diffx:difft')
      call getmem4d(dx%diffqx,jci1,jci2,ici1,ici2,1,kz,1,nqx,'diffx:diffqx')
    end subroutine allocate_diffx
!
    subroutine allocate_mod_atm_interface(ibltyp)
!
      implicit none
      integer(ik4) , intent(in) :: ibltyp
!
      call allocate_domain(mddom,.true.)

      call allocate_atmstate(atm1,ibltyp,.true.,one_exchange_point)
      call allocate_atmstate(atm2,ibltyp,.true.,one_exchange_point)
      call allocate_atmstate(atmx,ibltyp,.true.,one_exchange_point)
      call allocate_atmstate(atmc,ibltyp,.true.,zero_exchange_point)
      call allocate_atmstate(aten,ibltyp,.true.,zero_exchange_point)
      if ( ibltyp == 99 ) then
        call allocate_atmstate(holtten,ibltyp,.true.,zero_exchange_point)
      end if
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call allocate_atmstate(uwten,ibltyp,.true.,one_exchange_point)
      end if

      call allocate_surfstate(sfs,.true.)

      call allocate_slice(atms,ibltyp)

      call allocate_diffx(adf)

      call getmem2d(ts0,jce1,jce2,ice1,ice2,'atm_interface:ts0')
      call getmem2d(ts1,jce1,jce2,ice1,ice2,'atm_interface:ts1')

      call getmem3d(dstor,jde1,jde2,ide1,ide2,1,nsplit,'atm_interface:dstor')
      call getmem3d(hstor,jde1,jde2,ide1,ide2,1,nsplit,'atm_interface:hstor')
!
      call getmem2d(hgfact,jde1,jde2,ide1,ide2,'atm_interface:hgfact')
      call getmem2d(psdot,jde1,jde2,ide1,ide2,'atm_interface:psdot')
      call getmem3d(omega,jci1,jci2,ici1,ici2,1,kz,'atm_interface:omega')
      call getmem3d(qdot,jce1-ma%jbl1,jce2+ma%jbr1, &
                         ice1-ma%ibb1,ice2+ma%ibt1,1,kzp1,'atm_interface:qdot')

    end subroutine allocate_mod_atm_interface 
!
end module mod_atm_interface
