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
  use mod_constants , only : d_rfour
  use mod_runparams
  use mod_mppparam
  use mod_memutil
  use mpi

  private

  logical , public , parameter :: cross = .false.
  logical , public , parameter :: dot = .true.

!
! Storage for all the 3d prognostic variables in two
!     timesteps and all the 2d variables and constants
!
  type domain
    real(dp) , pointer , dimension(:,:) :: ht
    real(dp) , pointer , dimension(:,:) :: lndcat
    real(dp) , pointer , dimension(:,:) :: xlat
    real(dp) , pointer , dimension(:,:) :: xlon
    real(dp) , pointer , dimension(:,:) :: msfx
    real(dp) , pointer , dimension(:,:) :: msfd
    real(dp) , pointer , dimension(:,:) :: coriol
  end type domain

  type atmstate
    real(dp) , pointer , dimension(:,:,:) :: u
    real(dp) , pointer , dimension(:,:,:) :: v
    real(dp) , pointer , dimension(:,:,:) :: t
    real(dp) , pointer , dimension(:,:,:) :: qv
    real(dp) , pointer , dimension(:,:,:) :: qc
    real(dp) , pointer , dimension(:,:,:) :: tke
  end type atmstate

  type surfstate
    real(dp) , pointer , dimension(:,:) :: psa
    real(dp) , pointer , dimension(:,:) :: psb
    real(dp) , pointer , dimension(:,:) :: tga
    real(dp) , pointer , dimension(:,:) :: tgb
    real(dp) , pointer , dimension(:,:) :: rainc
    real(dp) , pointer , dimension(:,:) :: rainnc
    real(dp) , pointer , dimension(:,:) :: hfx
    real(dp) , pointer , dimension(:,:) :: qfx
    real(dp) , pointer , dimension(:,:) :: tgbb
    real(dp) , pointer , dimension(:,:) :: uvdrag
  end type surfstate

  type slice
    real(dp) , pointer , dimension(:,:,:) :: tb3d
    real(dp) , pointer , dimension(:,:,:) :: thx3d
    real(dp) , pointer , dimension(:,:,:) :: pb3d
    real(dp) , pointer , dimension(:,:,:) :: rhob3d
    real(dp) , pointer , dimension(:,:,:) :: ubx3d
    real(dp) , pointer , dimension(:,:,:) :: vbx3d
    real(dp) , pointer , dimension(:,:,:) :: ubd3d
    real(dp) , pointer , dimension(:,:,:) :: vbd3d
    real(dp) , pointer , dimension(:,:,:) :: rhb3d
    real(dp) , pointer , dimension(:,:,:) :: qvb3d
    real(dp) , pointer , dimension(:,:,:) :: qsb3d
    real(dp) , pointer , dimension(:,:,:) :: qcb3d
    real(dp) , pointer , dimension(:,:,:) :: tkeb3d
    real(dp) , pointer , dimension(:,:,:,:) :: chib3d
  end type slice

  type diffx
    real(dp) , pointer , dimension(:,:,:) :: difft
    real(dp) , pointer , dimension(:,:,:) :: difuu
    real(dp) , pointer , dimension(:,:,:) :: difuv
    real(dp) , pointer , dimension(:,:,:) :: diffq
  end type diffx

  type v3dbound
    real(dp) , pointer , dimension(:,:,:) :: b0
    real(dp) , pointer , dimension(:,:,:) :: b1
    real(dp) , pointer , dimension(:,:,:) :: bt
  end type v3dbound

  type v2dbound
    real(dp) , pointer , dimension(:,:) :: b0
    real(dp) , pointer , dimension(:,:) :: b1
    real(dp) , pointer , dimension(:,:) :: bt
  end type v2dbound

  type bound_area
    logical :: dotflag
    logical :: havebound
    logical , pointer , dimension(:,:) :: bsouth
    logical , pointer , dimension(:,:) :: bnorth
    logical , pointer , dimension(:,:) :: beast
    logical , pointer , dimension(:,:) :: bwest
    integer :: ns , nn , ne , nw
    integer :: nsp
    integer , pointer , dimension(:,:) :: ibnd
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
  public :: allocate_surfstate
  public :: deco1_bound , deco1_allocate_v3dbound , deco1_allocate_v2dbound
  public :: deco1_model

  real(dp) , public , pointer , dimension(:,:) :: hgfact
  real(dp) , public , pointer , dimension(:,:) :: psdot
  real(dp) , public , pointer, dimension(:,:,:) :: dstor
  real(dp) , public , pointer, dimension(:,:,:) :: hstor
!
  real(dp) , public , pointer , dimension(:,:,:) :: qdot , omega
!
  contains 
!
    subroutine deco1_model(lband)
      implicit none
      logical , intent(in) :: lband
      ma%bandflag = lband
      ma%hasleft  = .false.
      ma%hasright = .false.
      jde1  = 1
      jdi1  = 1
      jdii1 = 1
      jde2  = jxp
      jdi2  = jxp
      jdii2 = jxp
      ide1  = idot1
      idi1  = idot1+1
      idii1 = idot1+2
      ide2  = idot2
      idi2  = idot2-1
      idii2 = idot2-2
      if ( .not. lband ) then
        if ( myid == 0 ) then
          jde1  = jcross1
          jdi1  = jcross1+1
          jdii1 = jcross1+2
          ma%hasleft = .true.
        end if
        if ( myid == nproc-1 ) then
          jde2  = jxp
          jdi2  = jxp-1
          jdii2 = jxp-2
          ma%hasright = .true.
        end if
      end if
      jce1  = 1
      jci1  = 1
      jcii1 = 1
      jce2  = jxp
      jci2  = jxp
      jcii2 = jxp
      ice1  = icross1
      ici1  = icross1+1
      icii1 = icross1+2
      ice2  = icross2
      ici2  = icross2-1
      icii2 = icross2-2
      if ( .not. lband ) then
        if ( myid == 0 ) then
          jce1  = jcross1
          jci1  = jcross1+1
          jcii1 = jcross1+2
        end if
        if ( myid == nproc-1 ) then
          jce2  = jxp-1
          jci2  = jxp-2
          jcii2 = jxp-3
        end if
      end if
      ! In 1D deco, each processor ALWAYS HAS top and bottom.
      ma%hastop    = .true.
      ma%hasbottom = .true.
    end subroutine deco1_model

    subroutine deco1_bound(ldot,lband,ba)
      implicit none
      logical , intent(in) :: ldot , lband
      type(bound_area) , intent(out) :: ba
      integer :: ic
      integer :: igbb1 , igbb2 , igbt1 , igbt2
      integer :: jgbl1 , jgbl2 , jgbr1 , jgbr2
      integer :: i , j , i1 , i2 , j1 , j2 , iglob , jglob

      ba%dotflag = ldot
      call getmem2d(ba%ibnd,1,jxp,idot1,idot2,'deco1_bound:ibnd')
      call getmem2d(ba%bsouth,1,jxp,idot1,idot2,'deco1_bound:bsouth')
      call getmem2d(ba%bnorth,1,jxp,idot1,idot2,'deco1_bound:bnorth')
      call getmem2d(ba%beast,1,jxp,idot1,idot2,'deco1_bound:beast')
      call getmem2d(ba%bwest,1,jxp,idot1,idot2,'deco1_bound:bwest')
      if (ldot ) then
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
      i2 = iy
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
      if ( lband ) then
        ! Check for South boundary
        do i = i1 , i2
          iglob = i
          if ( iglob >= igbb1 .and. iglob <= igbb2 ) then
            do j = j1 , j2
              jglob = jxp*myid+j
              if ( jglob < jgbl1 .and. jglob > jgbr2 ) cycle
              ba%ibnd(j,i) = iglob-igbb1+2
              ba%bsouth(j,i) = .true.
              ba%ns = ba%ns+1
            end do
          end if
        end do
        ! North Boundary
        do i = i1 , i2
          iglob = i
          if ( iglob >= igbt1 .and. iglob <= igbt2 ) then
            do j = j1 , j2
              jglob = jxp*myid+j
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
          iglob = i
          if ( iglob >= igbb1 .and. iglob <= igbb2 ) then
            do j = j1 , j2
              jglob = jxp*myid+j
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
          iglob = i
          if ( iglob >= igbt1 .and. iglob <= igbt2 ) then
            do j = j1 , j2
              jglob = jxp*myid+j
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
          iglob = i
          if ( iglob < igbb1 .or. iglob > igbt2 ) cycle
          do j = j1 , j2
            jglob = jxp*myid+j
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
          iglob = i
          if ( iglob < igbb1 .or. iglob > igbt2 ) cycle
          do j = j1 , j2
            jglob = jxp*myid+j
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
    end subroutine deco1_bound

    subroutine deco1_allocate_v3dbound(xb,ke,ldot)
      implicit none
      type(v3dbound) , intent(out) :: xb
      integer , intent(in) :: ke
      logical , intent(in) :: ldot
      if ( ldot ) then
        call getmem3d(xb%b0,0,jxp+1,idot1,idot2,1,ke,'v3dbound:b0')
        call getmem3d(xb%b1,0,jxp+1,idot1,idot2,1,ke,'v3dbound:b1')
        call getmem3d(xb%bt,0,jxp+1,idot1,idot2,1,ke,'v3dbound:bt')
      else
        call getmem3d(xb%b0,0,jxp+1,icross1,icross2,1,ke,'v3dbound:b0')
        call getmem3d(xb%b1,0,jxp+1,icross1,icross2,1,ke,'v3dbound:b1')
        call getmem3d(xb%bt,0,jxp+1,icross1,icross2,1,ke,'v3dbound:bt')
      end if
    end subroutine deco1_allocate_v3dbound
!
    subroutine deco1_allocate_v2dbound(xb,ldot)
      implicit none
      type(v2dbound) , intent(out) :: xb
      logical , intent(in) :: ldot
      if ( ldot ) then
        call getmem2d(xb%b0,0,jxp+1,idot1,idot2,'v2dbound:b0')
        call getmem2d(xb%b1,0,jxp+1,idot1,idot2,'v2dbound:b1')
        call getmem2d(xb%bt,0,jxp+1,idot1,idot2,'v2dbound:bt')
      else
        call getmem2d(xb%b0,0,jxp+1,icross1,icross2,'v2dbound:b0')
        call getmem2d(xb%b1,0,jxp+1,icross1,icross2,'v2dbound:b1')
        call getmem2d(xb%bt,0,jxp+1,icross1,icross2,'v2dbound:bt')
      end if
    end subroutine deco1_allocate_v2dbound
!
    subroutine allocate_atmstate(atm,ibltyp,lpar,ib,jb)
      implicit none
      logical , intent(in) :: lpar
      integer , intent(in) :: ibltyp
      integer , intent(in) :: ib , jb
      type(atmstate) , intent(out) :: atm
      integer :: is , ie , js , je
      if (lpar) then
        is = idot1-ib
        ie = idot2+ib
        js = 1-jb
        je = jxp+jb
        call getmem3d(atm%u,js,je,is,ie,1,kz,'atmstate:u')
        call getmem3d(atm%v,js,je,is,ie,1,kz,'atmstate:v')
        call getmem3d(atm%t,js,je,is,ie,1,kz,'atmstate:t')
        call getmem3d(atm%qv,js,je,is,ie,1,kz,'atmstate:qv')
        call getmem3d(atm%qc,js,je,is,ie,1,kz,'atmstate:qc')
        if ( ibltyp == 2 .or. ibltyp == 99 ) then
          call getmem3d(atm%tke,js,je,is,ie,1,kzp1,'atmstate:tke')
        end if
      else
        call getmem3d(atm%u,jdot1,jdot2,idot1,idot2,1,kz,'atmstate:u')
        call getmem3d(atm%v,jdot1,jdot2,idot1,idot2,1,kz,'atmstate:v')
        call getmem3d(atm%t,jcross1,jcross2,icross1,icross2,1,kz,'atmstate:t')
        call getmem3d(atm%qv,jcross1,jcross2,icross1,icross2,1,kz,'atmstate:qv')
        call getmem3d(atm%qc,jcross1,jcross2,icross1,icross2,1,kz,'atmstate:qc')
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
        call getmem2d(dom%ht,0,jxp+1,idot1,idot2,'mod_atm_interface:ht')
        call getmem2d(dom%lndcat,1,jxp,idot1,idot2,'mod_atm_interface:lndcat')
        call getmem2d(dom%xlat,1,jxp,idot1,idot2,'mod_atm_interface:xlat')
        call getmem2d(dom%xlon,1,jxp,idot1,idot2,'mod_atm_interface:xlon')
        call getmem2d(dom%msfx,-1,jxp+2,idot1,idot2,'mod_atm_interface:msfx')
        call getmem2d(dom%msfd,-1,jxp+2,idot1,idot2,'mod_atm_interface:msfd')
        call getmem2d(dom%coriol,1,jxp,idot1,idot2,'mod_atm_interface:f')
      else
        call getmem2d(dom%ht,jdot1,jdot2,idot1,idot2,'atm_interface:ht')
        call getmem2d(dom%lndcat,jdot1,jdot2,idot1,idot2,'atm_interface:lndcat')
        call getmem2d(dom%xlat,jdot1,jdot2,idot1,idot2,'atm_interface:xlat')
        call getmem2d(dom%xlon,jdot1,jdot2,idot1,idot2,'atm_interface:xlon')
        call getmem2d(dom%msfx,jdot1,jdot2,idot1,idot2,'atm_interface:msfx')
        call getmem2d(dom%msfd,jdot1,jdot2,idot1,idot2,'atm_interface:msfd')
        call getmem2d(dom%coriol,jdot1,jdot2,idot1,idot2,'atm_interface:f')
      end if
    end subroutine allocate_domain
!
    subroutine allocate_surfstate(sfs,lpar)
      implicit none
      type(surfstate) , intent(out) :: sfs
      logical , intent(in) :: lpar
      if (lpar) then
        call getmem2d(sfs%psa,0,jxp+1,icross1,icross2,'surf:psa')
        call getmem2d(sfs%psb,0,jxp+1,icross1,icross2,'surf:psb')
        call getmem2d(sfs%tga,1,jxp+1,icross1,icross2,'surf:tga')
        call getmem2d(sfs%tgb,1,jxp+1,icross1,icross2,'surf:tgb')
        call getmem2d(sfs%hfx,1,jxp,icross1,icross2,'surf:hfx')
        call getmem2d(sfs%qfx,1,jxp,icross1,icross2,'surf:qfx')
        call getmem2d(sfs%rainc,1,jxp,icross1,icross2,'surf:rainc')
        call getmem2d(sfs%rainnc,1,jxp,icross1,icross2,'surf:rainnc')
        call getmem2d(sfs%tgbb,1,jxp,icross1,icross2,'surf:tgbb')
        call getmem2d(sfs%uvdrag,0,jxp,icross1,icross2,'surf:uvdrag')
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
      integer , intent(in) :: ibltyp
      call getmem3d(ax%pb3d,1,jxp,icross1,icross2,1,kz,'slice:pb3d')
      call getmem3d(ax%qsb3d,1,jxp,icross1,icross2,1,kz,'slice:qsb3d')
      call getmem3d(ax%rhb3d,1,jxp,icross1,icross2,1,kz,'slice:rhb3d')
      call getmem3d(ax%rhob3d,1,jxp,icross1,icross2,1,kz,'slice:rhob3d')
      call getmem3d(ax%ubx3d,-1,jxp+2,icross1,icross2,1,kz,'slice:ubx3d')
      call getmem3d(ax%vbx3d,-1,jxp+2,icross1,icross2,1,kz,'slice:vbx3d')
      call getmem3d(ax%thx3d,1,jxp,icross1,icross2,1,kz,'slice:thx3d')
      call getmem3d(ax%qcb3d,-1,jxp+2,icross1,icross2,1,kz,'slice:qcb3d')
      call getmem3d(ax%qvb3d,-1,jxp+2,icross1,icross2,1,kz,'slice:qvb3d')
      call getmem3d(ax%tb3d,-1,jxp+2,icross1,icross2,1,kz,'slice:tb3d')
      call getmem3d(ax%ubd3d,-1,jxp+2,idot1,idot2,1,kz,'slice:ubd3d')
      call getmem3d(ax%vbd3d,-1,jxp+2,idot1,idot2,1,kz,'slice:vbd3d')
      if ( ichem == 1 ) then
        call getmem4d(ax%chib3d,-1,jxp+2,icross1,icross2, &
                      1,kz,1,ntr,'slice:chib3d')
      end if
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call getmem3d(ax%tkeb3d,-1,jxp+2,icross1,icross2,1,kzp1,'slice:tkeb3d')
      end if
    end subroutine allocate_slice
!
    subroutine allocate_diffx(dx)
      implicit none
      type(diffx) , intent(out) :: dx
      call getmem3d(dx%difft,1,jxp,icross1,icross2,1,kz,'diffx:difft')
      call getmem3d(dx%difuu,1,jxp,idot1,idot2,1,kz,'diffx:difuu')
      call getmem3d(dx%difuv,1,jxp,idot1,idot2,1,kz,'diffx:difuv')
      call getmem3d(dx%diffq,1,jxp,icross1,icross2,1,kz,'diffx:diffq')
    end subroutine allocate_diffx
!
    subroutine allocate_mod_atm_interface(ibltyp)
!
      implicit none
      integer , intent(in) :: ibltyp
!
      call allocate_domain(mddom,.true.)

      call allocate_atmstate(atm1,ibltyp,.true.,0,2)
      call allocate_atmstate(atm2,ibltyp,.true.,0,2)
      call allocate_atmstate(atmx,ibltyp,.true.,0,1)
      call allocate_atmstate(atmc,ibltyp,.true.,0,0)
      call allocate_atmstate(aten,ibltyp,.true.,0,0)
      if ( ibltyp == 99 ) then
        call allocate_atmstate(holtten,ibltyp,.true.,0,0)
      end if
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        call allocate_atmstate(uwten,ibltyp,.true.,0,1)
      end if

      call allocate_surfstate(sfs,.true.)

      call allocate_slice(atms,ibltyp)

      call allocate_diffx(adf)

      call getmem3d(dstor,1,jxp,idot1,idot2,1,nsplit,'mod_atm_interface:dstor')
      call getmem3d(hstor,1,jxp,idot1,idot2,1,nsplit,'mod_atm_interface:hstor')
!
      call getmem2d(hgfact,1,jxp,idot1,idot2,'mod_atm_interface:hgfact')
      call getmem2d(psdot,1,jxp,idot1,idot2,'mod_atm_interface:psdot')
      call getmem3d(omega,1,jxp,icross1,icross2,1,kz,'mod_atm_interface:omega')
      call getmem3d(qdot,0,jxp+1,idot1,idot2,1,kzp1,'mod_atm_interface:qdot')

    end subroutine allocate_mod_atm_interface 
!
end module mod_atm_interface
