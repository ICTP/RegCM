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

module mod_atm_stub

  use mod_dynparam
  use mod_stdio
  use mod_date
  use mod_constants
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_service
  use mod_memutil
  use mod_regcm_types
  use mod_stdatm
  use mod_zita
  use mod_ncio

  implicit none (type, external)

  private

  type(domain), public :: mddom
  type(domain_subgrid), public :: mdsub
  type(surfstate), public :: sfs

  type(v3dbound), public :: xtb, xqb, xub, xvb, xlb, xib, xppb, xwwb
  type(v2dbound), public :: xpsb, xtsb, xpaib

  type(lm_exchange), public :: lm
  type(lm_state), public :: lms

  integer(ik4) :: ix1, ix2, jx1, jx2
  integer(ik4) :: id1, id2, jd1, jd2

  real(rkx), public :: rdnnsg
  real(rkx) :: rdtbdy

  logical, parameter :: cross = .false.
  logical, parameter :: dot = .true.

  real(rkx), dimension(:,:), pointer, contiguous, public :: patm
  real(rkx), dimension(:,:), pointer, contiguous, public :: tatm
  real(rkx), dimension(:,:), pointer, contiguous, public :: uatm
  real(rkx), dimension(:,:), pointer, contiguous, public :: vatm
  real(rkx), dimension(:,:), pointer, contiguous, public :: thatm
  real(rkx), dimension(:,:), pointer, contiguous, public :: qvatm
  real(rkx), dimension(:,:), pointer, contiguous, public :: zatm
  real(rkx), dimension(:,:), pointer, contiguous, public :: rho
  real(rkx), dimension(:,:), pointer, contiguous, public :: ps
  real(rkx), dimension(:,:), pointer, contiguous, public :: tp
  real(rkx), dimension(:,:), pointer, contiguous, public :: coszrs
  real(rkx), dimension(:,:), pointer, contiguous, public :: fsw
  real(rkx), dimension(:,:), pointer, contiguous, public :: flw
  real(rkx), dimension(:,:), pointer, contiguous, public :: flwd
  real(rkx), dimension(:,:), pointer, contiguous, public :: solar
  real(rkx), dimension(:,:), pointer, contiguous, public :: pptc
  real(rkx), dimension(:,:), pointer, contiguous, public :: pptnc
  real(rkx), dimension(:,:), pointer, contiguous, public :: swdir
  real(rkx), dimension(:,:), pointer, contiguous, public :: swdif
  real(rkx), dimension(:,:), pointer, contiguous, public :: lwdir
  real(rkx), dimension(:,:), pointer, contiguous, public :: lwdif
  real(rkx), dimension(:,:), pointer, contiguous, public :: totc
  real(rkx), dimension(:,:), pointer, contiguous, public :: zeta, fmzf

  public :: allocate_mod_atm_interface, allocate_surface_model
  public :: setup_model_indexes
  public :: init_bdy, bdyin, atmval

  interface timeint
    module procedure timeint2, timeint3
  end interface timeint

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
        call getmem(xb%b0,jce1,jce2,ice1,ice2,1,ke,'v3dbound:b0')
        call getmem(xb%b1,jce1,jce2,ice1,ice2,1,ke,'v3dbound:b1')
        call getmem(xb%bt,jce1,jce2,ice1,ice2,1,ke,'v3dbound:bt')
      end if
    end subroutine allocate_v3dbound

    subroutine allocate_v2dbound(xb,ldot)
      implicit none (type, external)
      type(v2dbound), intent(inout) :: xb
      logical, intent(in) :: ldot
      if ( ldot ) then
        call getmem(xb%b0,jde1,jde2,ide1,ide2,'v2dbound:b0')
        call getmem(xb%b1,jde1,jde2,ide1,ide2,'v2dbound:b1')
        call getmem(xb%bt,jde1,jde2,ide1,ide2,'v2dbound:bt')
      else
        call getmem(xb%b0,jce1,jce2,ice1,ice2,'v2dbound:b0')
        call getmem(xb%b1,jce1,jce2,ice1,ice2,'v2dbound:b1')
        call getmem(xb%bt,jce1,jce2,ice1,ice2,'v2dbound:bt')
      end if
    end subroutine allocate_v2dbound

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
      call getmem(sfs%psa,jce1,jce2,ice1,ice2,'surf:psa')
      call getmem(sfs%psdota,jde1,jde2,ide1,ide2,'surf:psdota')
      call getmem(sfs%psb,jx1,jx2,ix1,ix2,'surf:psb')
      call getmem(sfs%psdotb,jd1,jd2,id1,id2,'surf:psdotb')
      call getmem(sfs%psc,jce1,jce2,ice1,ice2,'surf:psc')
      call getmem(sfs%tg,jci1,jci2,ici1,ici2,'surf:tg')
      call getmem(sfs%hfx,jci1,jci2,ici1,ici2,'surf:hfx')
      call getmem(sfs%qfx,jci1,jci2,ici1,ici2,'surf:qfx')
      call getmem(sfs%rainc,jci1,jci2,ici1,ici2,'surf:rainc')
      call getmem(sfs%rainnc,jci1,jci2,ici1,ici2,'surf:rainnc')
      call getmem(sfs%tgbb,jci1,jci2,ici1,ici2,'surf:tgbb')
      call getmem(sfs%uvdrag,jci1,jci2,ici1,ici2,'surf:uvdrag')
      call getmem(sfs%zo,jci1,jci2,ici1,ici2,'surf:zo')
      call getmem(sfs%ram1,jci1,jci2,ici1,ici2,'surf:ram1')
      call getmem(sfs%rah1,jci1,jci2,ici1,ici2,'surf:rah1')
      call getmem(sfs%br,jci1,jci2,ici1,ici2,'surf:br')
      call getmem(sfs%q2m,jci1,jci2,ici1,ici2,'surf:q2m')
      call getmem(sfs%ustar,jci1,jci2,ici1,ici2,'surf:ustar')
      call getmem(sfs%w10m,jci1,jci2,ici1,ici2,'surf:w10m')
      call getmem(sfs%u10m,jci1,jci2,ici1,ici2,'surf:u10m')
      call getmem(sfs%v10m,jci1,jci2,ici1,ici2,'surf:v10m')
    end subroutine allocate_surfstate

    subroutine allocate_mod_atm_interface
      implicit none (type, external)

      id1 = ide1ga
      id2 = ide2ga
      jd1 = jde1ga
      jd2 = jde2ga
      ix1 = ice1ga
      ix2 = ice2ga
      jx1 = jce1ga
      jx2 = jce2ga

      call allocate_domain(mddom)
      call allocate_domain_subgrid(mdsub)
      call allocate_surfstate(sfs)
      call allocate_v2dbound(xpsb,cross)
      call allocate_v2dbound(xtsb,cross)
      call allocate_v3dbound(xtb,kz,cross)
      call allocate_v3dbound(xqb,kz,cross)
      call allocate_v3dbound(xub,kz,dot)
      call allocate_v3dbound(xlb,kz,dot)
      call allocate_v3dbound(xib,kz,dot)
      call allocate_v3dbound(xvb,kz,dot)
      if ( idynamic == 3 ) then
        call allocate_v2dbound(xpaib,cross)
      end if
      if ( idynamic == 3 ) then
        call getmem(zeta,jce1,jce2,ice1,ice2,'lm:zeta')
        call getmem(fmzf,jce1,jce2,ice1,ice2,'lm:fmzf')
      end if
    end subroutine allocate_mod_atm_interface

    subroutine allocate_surface_model
      implicit none (type, external)

      rdnnsg = d_one/real(nnsg,rkx)

      call getmem(lms%sent,1,nnsg,jci1,jci2,ici1,ici2,'lm:sent')
      call getmem(lms%evpr,1,nnsg,jci1,jci2,ici1,ici2,'lm:evpr')
      call getmem(lms%deltat,1,nnsg,jci1,jci2,ici1,ici2,'lm:deltat')
      call getmem(lms%deltaq,1,nnsg,jci1,jci2,ici1,ici2,'lm:deltaq')
      call getmem(lms%drag,1,nnsg,jci1,jci2,ici1,ici2,'lm:drag')
      call getmem(lms%ustar,1,nnsg,jci1,jci2,ici1,ici2,'lm:ustar')
      call getmem(lms%w10m,1,nnsg,jci1,jci2,ici1,ici2,'lm:w10m')
      call getmem(lms%zo,1,nnsg,jci1,jci2,ici1,ici2,'lm:zo')
      call getmem(lms%rhoa,1,nnsg,jci1,jci2,ici1,ici2,'lm:rho')
      call getmem(lms%lncl,1,nnsg,jci1,jci2,ici1,ici2,'lm:lncl')
      call getmem(lms%prcp,1,nnsg,jci1,jci2,ici1,ici2,'lm:prcp')
      call getmem(lms%snwm,1,nnsg,jci1,jci2,ici1,ici2,'lm:snwm')
      call getmem(lms%trnof,1,nnsg,jci1,jci2,ici1,ici2,'lm:trnof')
      call getmem(lms%srnof,1,nnsg,jci1,jci2,ici1,ici2,'lm:srnof')
      call getmem(lms%xlai,1,nnsg,jci1,jci2,ici1,ici2,'lm:xlai')
      call getmem(lms%sfcp,1,nnsg,jci1,jci2,ici1,ici2,'lm:sfcp')
      call getmem(lms%q2m,1,nnsg,jci1,jci2,ici1,ici2,'lm:q2m')
      call getmem(lms%t2m,1,nnsg,jci1,jci2,ici1,ici2,'lm:t2m')
      call getmem(lms%u10m,1,nnsg,jci1,jci2,ici1,ici2,'lm:u10m')
      call getmem(lms%v10m,1,nnsg,jci1,jci2,ici1,ici2,'lm:v10m')
      call getmem(lms%ram1,1,nnsg,jci1,jci2,ici1,ici2,'lm:ram1')
      call getmem(lms%rah1,1,nnsg,jci1,jci2,ici1,ici2,'lm:rah1')
      call getmem(lms%br,1,nnsg,jci1,jci2,ici1,ici2,'lm:br')
      call getmem(lms%taux,1,nnsg,jci1,jci2,ici1,ici2,'lm:taux')
      call getmem(lms%tauy,1,nnsg,jci1,jci2,ici1,ici2,'lm:tauy')
      call getmem(lms%wt,1,nnsg,jci1,jci2,ici1,ici2,'lm:wt')
      call getmem(lms%swalb,1,nnsg,jci1,jci2,ici1,ici2,'lm:swalb')
      call getmem(lms%lwalb,1,nnsg,jci1,jci2,ici1,ici2,'lm:lwalb')
      call getmem(lms%swdiralb,1,nnsg,jci1,jci2,ici1,ici2,'lm:swdiralb')
      call getmem(lms%lwdiralb,1,nnsg,jci1,jci2,ici1,ici2,'lm:lwdiralb')
      call getmem(lms%swdifalb,1,nnsg,jci1,jci2,ici1,ici2,'lm:swdifalb')
      call getmem(lms%lwdifalb,1,nnsg,jci1,jci2,ici1,ici2,'lm:lwdifalb')

      call getmem(lms%gwet,1,nnsg,jci1,jci2,ici1,ici2,'lm:gwet')
      call getmem(lms%sw,1,nnsg,jci1,jci2,ici1,ici2,1,num_soil_layers,'lm:sw')

      call assignpnt(lms%sw,lms%ssw,1)
      call assignpnt(lms%sw,lms%rsw,2)
      call assignpnt(lms%sw,lms%tsw,3)
      call getmem(lms%tgbb,1,nnsg,jci1,jci2,ici1,ici2,'lm:tgbb')
      call getmem(lms%tgrd,1,nnsg,jci1,jci2,ici1,ici2,'lm:tgrd')
      call getmem(lms%tgbrd,1,nnsg,jci1,jci2,ici1,ici2,'lm:tgbrd')
      call getmem(lms%tlef,1,nnsg,jci1,jci2,ici1,ici2,'lm:tlef')
      call getmem(lms%taf,1,nnsg,jci1,jci2,ici1,ici2,'lm:taf')
      call getmem(lms%sigf,1,nnsg,jci1,jci2,ici1,ici2,'lm:sigf')
      call getmem(lms%sfice,1,nnsg,jci1,jci2,ici1,ici2,'lm:sfice')
      call getmem(lms%snag,1,nnsg,jci1,jci2,ici1,ici2,'lm:snag')
      call getmem(lms%ldew,1,nnsg,jci1,jci2,ici1,ici2,'lm:ldew')
      call getmem(lms%sncv,1,nnsg,jci1,jci2,ici1,ici2,'lm:sncv')
      call getmem(lms%scvk,1,nnsg,jci1,jci2,ici1,ici2,'lm:scvk')
      call getmem(lms%um10,1,nnsg,jci1,jci2,ici1,ici2,'lm:um10')
      call getmem(lms%emisv,1,nnsg,jci1,jci2,ici1,ici2,'lm:emisv')
      call getmem(lms%vocemiss,1,nnsg,jci1,jci2,ici1,ici2,1,ntr,'lm:vocemiss')
      call getmem(lms%dustemiss,1,nnsg,jci1,jci2,ici1,ici2,1,4,'lm:dustemiss')
      call getmem(lms%ddepv,1,nnsg,jci1,jci2,ici1,ici2,1,ntr,'lm:ddepv')
      call getmem(lms%sw_vol,1,nnsg,jci1,jci2, &
                                      ici1,ici2,1,num_soil_layers,'lm:sw_vol')
      call getmem(lms%tsoi,1,nnsg,jci1,jci2, &
                                    ici1,ici2,1,num_soil_layers,'lm:tsoi')

      call getmem(tatm,jci1,jci2,ici1,ici2,'lm:tatm')
      call getmem(patm,jci1,jci2,ici1,ici2,'lm:patm')
      call getmem(uatm,jci1,jci2,ici1,ici2,'lm:uatm')
      call getmem(vatm,jci1,jci2,ici1,ici2,'lm:vatm')
      call getmem(thatm,jci1,jci2,ici1,ici2,'lm:thatm')
      call getmem(qvatm,jci1,jci2,ici1,ici2,'lm:qvatm')
      call getmem(zatm,jci1,jci2,ici1,ici2,'lm:zatm')
      call getmem(rho,jci1,jci2,ici1,ici2,'lm:rho')
      call getmem(ps,jci1,jci2,ici1,ici2,'lm:ps')
      call getmem(tp,jci1,jci2,ici1,ici2,'lm:tp')
      call getmem(coszrs,jci1,jci2,ici1,ici2,'lm:coszrs')
      call getmem(flw,jci1,jci2,ici1,ici2,'lm:flw')
      call getmem(fsw,jci1,jci2,ici1,ici2,'lm:fsw')
      call getmem(flwd,jci1,jci2,ici1,ici2,'lm:flwd')
      call getmem(solar,jci1,jci2,ici1,ici2,'lm:solar')
      call getmem(pptc,jci1,jci2,ici1,ici2,'lm:pptc')
      call getmem(pptnc,jci1,jci2,ici1,ici2,'lm:pptnc')
      call getmem(swdir,jci1,jci2,ici1,ici2,'lm:swdir')
      call getmem(swdif,jci1,jci2,ici1,ici2,'lm:swdif')
      call getmem(lwdir,jci1,jci2,ici1,ici2,'lm:lwdir')
      call getmem(lwdif,jci1,jci2,ici1,ici2,'lm:lwdif')
      call getmem(totc,jci1,jci2,ici1,ici2,'lm:totc')
      call cl_setup(lndcomm,mddom%mask,mdsub%mask)
      call cl_setup(ocncomm,mddom%mask,mdsub%mask,.true.)

#ifdef DEBUG
      write(ndebug,*) 'TOTAL POINTS FOR LAND  IN LNDCOMM : ', &
        lndcomm%linear_npoint_sg(myid+1)
      write(ndebug,*) 'Cartesian p ', lndcomm%cartesian_npoint_g
      write(ndebug,*) 'Cartesian d ', lndcomm%cartesian_displ_g
      write(ndebug,*) 'Linear    p ', lndcomm%linear_npoint_g
      write(ndebug,*) 'Linear    d ', lndcomm%linear_displ_g
      write(ndebug,*) 'Subgrid Cartesian p ', lndcomm%cartesian_npoint_sg
      write(ndebug,*) 'Subgrid Cartesian d ', lndcomm%cartesian_displ_sg
      write(ndebug,*) 'Subgrid Linear    p ', lndcomm%linear_npoint_sg
      write(ndebug,*) 'Subgrid Linear    d ', lndcomm%linear_displ_sg
      write(ndebug,*) 'TOTAL POINTS FOR OCEAN IN OCNCOMM : ', &
        ocncomm%linear_npoint_sg(myid+1)
      write(ndebug,*) 'Cartesian p ', ocncomm%cartesian_npoint_g
      write(ndebug,*) 'Cartesian d ', ocncomm%cartesian_displ_g
      write(ndebug,*) 'Linear    p ', ocncomm%linear_npoint_g
      write(ndebug,*) 'Linear    d ', ocncomm%linear_displ_g
      write(ndebug,*) 'Subgrid Cartesian p ', ocncomm%cartesian_npoint_sg
      write(ndebug,*) 'Subgrid Cartesian d ', ocncomm%cartesian_displ_sg
      write(ndebug,*) 'Subgrid Linear    p ', ocncomm%linear_npoint_sg
      write(ndebug,*) 'Subgrid Linear    d ', ocncomm%linear_displ_sg
#endif

      call assignpnt(mddom%xlat,lm%xlat)
      call assignpnt(mddom%xlon,lm%xlon)
      call assignpnt(mddom%lndcat,lm%lndcat)
      call assignpnt(mddom%ldmsk,lm%ldmsk)
      call assignpnt(mddom%iveg,lm%iveg)
      call assignpnt(mddom%itex,lm%itex)
      call assignpnt(mddom%ht,lm%ht)
      call assignpnt(mddom%snowam,lm%snowam)
      call assignpnt(mddom%smoist,lm%smoist)
      call assignpnt(mddom%rmoist,lm%rmoist)
      call assignpnt(mddom%rts,lm%rts)
      call assignpnt(mdsub%xlat,lm%xlat1)
      call assignpnt(mdsub%xlon,lm%xlon1)
      call assignpnt(mdsub%area,lm%area1)
      call assignpnt(mdsub%lndcat,lm%lndcat1)
      call assignpnt(mdsub%ldmsk,lm%ldmsk1)
      call assignpnt(mdsub%ht,lm%ht1)
      call assignpnt(mdsub%iveg,lm%iveg1)
      call assignpnt(mdsub%itex,lm%itex1)
      call assignpnt(mdsub%dhlake,lm%dhlake1)
      call assignpnt(uatm,lm%uatm)
      call assignpnt(vatm,lm%vatm)
      call assignpnt(thatm,lm%thatm)
      call assignpnt(tatm,lm%tatm)
      call assignpnt(patm,lm%patm)
      call assignpnt(qvatm,lm%qvatm)
      call assignpnt(zatm,lm%hgt)

      call assignpnt(rho,lm%rhox)
      call assignpnt(ps,lm%sfps)
      call assignpnt(tp,lm%sfta)

      call assignpnt(sfs%hfx,lm%hfx)
      call assignpnt(sfs%qfx,lm%qfx)
      call assignpnt(sfs%uvdrag,lm%uvdrag)
      call assignpnt(sfs%ustar,lm%ustar)
      call assignpnt(sfs%w10m,lm%w10m)
      call assignpnt(sfs%zo,lm%zo)
      call assignpnt(sfs%tg,lm%tg)
      call assignpnt(sfs%tgbb,lm%tgbb)
      call assignpnt(sfs%u10m,lm%u10m)
      call assignpnt(sfs%v10m,lm%v10m)
      call assignpnt(sfs%ram1,lm%ram1)
      call assignpnt(sfs%rah1,lm%rah1)
      call assignpnt(sfs%br,lm%br)
      call assignpnt(sfs%q2m,lm%q2m)
      call assignpnt(pptc,lm%cprate)
      call assignpnt(pptnc,lm%ncprate)
      call assignpnt(coszrs,lm%zencos)
      call assignpnt(fsw,lm%rswf)
      call assignpnt(flw,lm%rlwf)
      call assignpnt(flwd,lm%dwrlwf)
      call assignpnt(solar,lm%solar)
      call assignpnt(swdir,lm%swdir)
      call assignpnt(swdif,lm%swdif)
      call assignpnt(lwdir,lm%lwdir)
      call assignpnt(lwdif,lm%lwdif)

    end subroutine allocate_surface_model

    subroutine init_bdy
      implicit none (type, external)
      character(len=32) :: appdat
      type (rcm_time_and_date) :: icbc_date
      integer(ik4) :: i, j, datefound
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'init_bdy'
      integer(ik4), save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif

      rdtbdy = d_one / dtbdys
      bdydate1 = idate1
      bdydate2 = idate1
      xbctime = d_zero

      if ( bdydate1 == globidate1 ) then
        icbc_date = bdydate1
      else
        icbc_date = monfirst(bdydate1)
      end if

      call open_icbc(icbc_date)
      if ( sfbcread == 1 ) then
        call open_clmbc(icbc_date)
      end if

      datefound = icbc_search(bdydate1)
      if (datefound < 0) then
        !
        ! Cannot run without initial conditions
        !
        appdat = tochar(bdydate1)
        call fatal(__FILE__,__LINE__,'ICBC for '//appdat//' not found')
      end if

      call read_icbc(xpsb%b0,xtsb%b0,mddom%ldmsk,xub%b0,xvb%b0,xtb%b0, &
                     xqb%b0,xlb%b0,xib%b0,xppb%b0,xwwb%b0)

      if ( myid == italk ) then
        appdat = tochar(bdydate1)
        if ( rcmtimer%start( ) ) then
          write(stdout,*) 'READY IC DATA for ', appdat
        else
          write(stdout,*) 'READY BC DATA for ', appdat
        end if
      end if

      if ( idynamic == 3 ) then
        xpsb%b0(:,:) = xpsb%b0(:,:)*d_100
        call exchange_lr(xub%b0,1,jde1,jde2,ide1,ide2,1,kz)
        call exchange_bt(xvb%b0,1,jde1,jde2,ide1,ide2,1,kz)
      else
        xpsb%b0(:,:) = (xpsb%b0(:,:)*d_r10)-ptop
        call exchange(xub%b0,1,jde1,jde2,ide1,ide2,1,kz)
        call exchange(xvb%b0,1,jde1,jde2,ide1,ide2,1,kz)
      end if

      bdydate2 = bdydate2 + intbdy
      if ( myid == italk ) then
        write(stdout,'(a,a,a,i8)') ' SEARCH BC data for ', tochar10(bdydate2), &
                        ', step = ', rcmtimer%lcount
      end if
      datefound = icbc_search(bdydate2)
      if ( datefound < 0 ) then
        call open_icbc(monfirst(bdydate2))
        datefound = icbc_search(bdydate2)
        if ( datefound < 0 ) then
          appdat = tochar(bdydate2)
          call fatal(__FILE__,__LINE__,'ICBC for '//appdat//' not found')
        end if
      end if

      call read_icbc(xpsb%b1,xtsb%b1,mddom%ldmsk,xub%b1,xvb%b1,xtb%b1, &
                     xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1)

      if ( myid == italk ) then
        write (stdout,*) 'READY  BC from     ', &
              tochar10(bdydate1), ' to ', tochar10(bdydate2)
      end if

      bdydate1 = bdydate2
      !
      ! Repeat for T2
      !
      if ( idynamic == 3 ) then
        xpsb%b1(:,:) = xpsb%b1(:,:)*d_100
        call exchange_lr(xub%b1,1,jde1,jde2,ide1,ide2,1,kz)
        call exchange_bt(xvb%b1,1,jde1,jde2,ide1,ide2,1,kz)
      else
        xpsb%b1(:,:) = (xpsb%b1(:,:)*d_r10)-ptop
        call exchange(xub%b1,1,jde1,jde2,ide1,ide2,1,kz)
        call exchange(xvb%b1,1,jde1,jde2,ide1,ide2,1,kz)
      end if
      !
      ! Calculate time varying component
      !
      call timeint(xub%b1,xub%b0,xub%bt,jde1,jde2,ide1,ide2,1,kz)
      call timeint(xvb%b1,xvb%b0,xvb%bt,jde1,jde2,ide1,ide2,1,kz)
      call timeint(xtb%b1,xtb%b0,xtb%bt,jce1,jce2,ice1,ice2,1,kz)
      call timeint(xqb%b1,xqb%b0,xqb%bt,jce1,jce2,ice1,ice2,1,kz)
      call timeint(xtsb%b1,xtsb%b0,xtsb%bt,jce1,jce2,ice1,ice2)
      call timeint(xpsb%b1,xpsb%b0,xpsb%bt,jce1,jce2,ice1,ice2)
      if ( idynamic == 3 ) then
        call paicompute(xpsb%b0,xtb%b0,xqb%b0,xpaib%b0)
        call paicompute(xpsb%b1,xtb%b1,xqb%b1,xpaib%b1)
        call timeint(xpaib%b1,xpaib%b0,xpaib%bt,jce1,jce2,ice1,ice2)
        call exchange_lr(xub%bt,1,jde1,jde2,ide1,ide2,1,kz)
        call exchange_bt(xvb%bt,1,jde1,jde2,ide1,ide2,1,kz)
      else
        call exchange(xub%bt,1,jde1,jde2,ide1,ide2,1,kz)
        call exchange(xvb%bt,1,jde1,jde2,ide1,ide2,1,kz)
      end if

      if ( rcmtimer%start( ) ) then
        do i = ici1, ici2
          do j = jci1, jci2
            sfs%tg(j,i) = xtsb%b0(j,i)
          end do
        end do
      end if

#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine init_bdy

    subroutine bdyin
      implicit none (type, external)
      integer(ik4) :: datefound
      character(len=32) :: appdat
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'bdyin'
      integer(ik4), save :: idindx = 0
      call time_begin(subroutine_name,idindx)
#endif

      xbctime = d_zero

      xub%b0(:,:,:) = xub%b1(:,:,:)
      xvb%b0(:,:,:) = xvb%b1(:,:,:)
      xtb%b0(:,:,:) = xtb%b1(:,:,:)
      xqb%b0(:,:,:) = xqb%b1(:,:,:)
      xtsb%b0(:,:) = xtsb%b1(:,:)
      xpsb%b0(:,:) = xpsb%b1(:,:)
      if ( idynamic == 3 ) then
        xpaib%b0(:,:) = xpaib%b1(:,:)
      end if

      bdydate2 = bdydate2 + intbdy
      if ( myid == italk ) then
        write(stdout,'(a,a,a,i8)') ' SEARCH BC data for ', tochar10(bdydate2), &
                        ', step = ', rcmtimer%lcount
      end if
      datefound = icbc_search(bdydate2)
      if ( datefound < 0 ) then
        call open_icbc(monfirst(bdydate2))
        datefound = icbc_search(bdydate2)
        if ( datefound < 0 ) then
          appdat = tochar(bdydate2)
          call fatal(__FILE__,__LINE__,'ICBC for '//appdat//' not found')
        end if
      end if
      call read_icbc(xpsb%b1,xtsb%b1,mddom%ldmsk,xub%b1,xvb%b1,xtb%b1, &
                     xqb%b1,xlb%b1,xib%b1,xppb%b1,xwwb%b1)

      if ( idynamic == 3 ) then
        xpsb%b1(:,:) = xpsb%b1(:,:)*d_100
        call exchange_lr(xub%b1,1,jde1,jde2,ide1,ide2,1,kz)
        call exchange_bt(xvb%b1,1,jde1,jde2,ide1,ide2,1,kz)
      else
        xpsb%b1(:,:) = (xpsb%b1(:,:)*d_r10)-ptop
        call exchange(xub%b1,1,jde1,jde2,ide1,ide2,1,kz)
        call exchange(xvb%b1,1,jde1,jde2,ide1,ide2,1,kz)
      end if
      call timeint(xpsb%b1,xpsb%b0,xpsb%bt,jce1,jce2,ice1,ice2)

      ! Linear time interpolation

      call timeint(xub%b1,xub%b0,xub%bt,jde1,jde2,ide1,ide2,1,kz)
      call timeint(xvb%b1,xvb%b0,xvb%bt,jde1,jde2,ide1,ide2,1,kz)
      call timeint(xtb%b1,xtb%b0,xtb%bt,jce1,jce2,ice1,ice2,1,kz)
      call timeint(xqb%b1,xqb%b0,xqb%bt,jce1,jce2,ice1,ice2,1,kz)
      call timeint(xtsb%b1,xtsb%b0,xtsb%bt,jce1,jce2,ice1,ice2)

      if ( idynamic == 3 ) then
        call paicompute(xpsb%b0,xtb%b0,xqb%b0,xpaib%b0)
        call paicompute(xpsb%b1,xtb%b1,xqb%b1,xpaib%b1)
        call timeint(xpaib%b1,xpaib%b0,xpaib%bt,jce1,jce2,ice1,ice2)
        call exchange_lr(xub%bt,1,jde1,jde2,ide1,ide2,1,kz)
        call exchange_bt(xvb%bt,1,jde1,jde2,ide1,ide2,1,kz)
      else
        call exchange(xub%bt,1,jde1,jde2,ide1,ide2,1,kz)
        call exchange(xvb%bt,1,jde1,jde2,ide1,ide2,1,kz)
      end if

      if ( myid == italk ) then
        write (stdout,*) 'READY  BC from     ', &
              tochar10(bdydate1), ' to ', tochar10(bdydate2)
      end if

      bdydate1 = bdydate2

#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine bdyin

    subroutine atmval
      implicit none (type, external)
      integer(ik4) :: i, j
      real(rkx) :: psb, cell, zq, xt

      xt = xbctime

      do i = ici1, ici2
        do j = jci1, jci2
          tatm(j,i) = xtb%b0(j,i,kz) + xt*xtb%bt(j,i,kz)
          if ( idynamic == 3 ) then
            ps(j,i) = xpsb%b0(j,i) + xt*xpsb%bt(j,i)
            psb = xpaib%b0(j,i) + xt*xpaib%bt(j,i)
            patm(j,i) = (psb**cpovr) * p00
            uatm(j,i) = 0.50_rkx*(xub%b0(j,i,kz) + xub%b0(j+1,i,kz) + &
                              xt*(xub%bt(j,i,kz) + xub%bt(j+1,i,kz)))
            vatm(j,i) = 0.50_rkx*(xvb%b0(j,i,kz) + xvb%b0(j,i+1,kz) + &
                              xt*(xvb%bt(j,i,kz) + xvb%bt(j,i+1,kz)))
          else
            psb = xpsb%b0(j,i) + xt*xpsb%bt(j,i)
            cell = ptop /psb
            ps(j,i) = (psb + ptop)*d_1000
            patm(j,i) = (hsigma(kz)*psb + ptop)*d_1000
            uatm(j,i) = 0.25_rkx*(xub%b0(j,i,kz) + xub%b0(j+1,i,kz) + &
                                  xub%b0(j,i+1,kz) + xub%b0(j+1,i+1,kz) + &
                              xt*(xub%bt(j,i,kz) + xub%bt(j+1,i,kz) + &
                                  xub%bt(j,i+1,kz) + xub%bt(j+1,i+1,kz)))
            vatm(j,i) = 0.25_rkx*(xvb%b0(j,i,kz) + xvb%b0(j+1,i,kz) + &
                                  xvb%b0(j,i+1,kz) + xvb%b0(j+1,i+1,kz) + &
                              xt*(xvb%bt(j,i,kz) + xvb%bt(j+1,i,kz) + &
                                  xvb%bt(j,i+1,kz) + xvb%bt(j+1,i+1,kz)))
          end if
          qvatm(j,i) = xqb%b0(j,i,kz) + xt*xqb%bt(j,i,kz)
          rho(j,i) = ps(j,i)/(rgas*tatm(j,i))
          thatm(j,i) = tatm(j,i) / (p00/patm(j,i))**rovcp
          tp(j,i) = tatm(j,i) / (ps(j,i)/patm(j,i))**rovcp
          if ( idynamic == 3 ) then
            zatm(j,i) = zeta(j,i)
          else
            zq = rovg * tatm(j,i) * log((sigma(kzp1)+cell)/(sigma(kz)+cell))
            zatm(j,i) = d_half*zq
          end if
        end do
      end do
      xbctime = xbctime + dtsrf
    end subroutine atmval

    subroutine timeint2(a,b,c,j1,j2,i1,i2)
      implicit none (type, external)
      real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: a, b
      real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: c
      integer(ik4), intent(in) :: j1, j2, i1, i2
      integer(ik4) :: i, j
      do i = i1, i2
        do j = j1, j2
          c(j,i) = (a(j,i)-b(j,i))*rdtbdy
        end do
      end do
    end subroutine timeint2

    subroutine timeint3(a,b,c,j1,j2,i1,i2,k1,k2)
      implicit none (type, external)
      real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: a, b
      real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: c
      integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
      integer(ik4) :: i, j, k
      do k = k1, k2
        do i = i1, i2
          do j = j1, j2
            c(j,i,k) = (a(j,i,k)-b(j,i,k))*rdtbdy
          end do
        end do
      end do
    end subroutine timeint3

  subroutine paicompute(xpsb,xtb,xqb,xpaib)
    implicit none (type, external)
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: xpsb
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: xtb, xqb
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: xpaib
    real(rkx) :: p, tv, zz1, zz2, zlr
    integer(ik4) :: i, j
    ! Hydrostatic initialization of pai
    do i = ice1, ice2
      do j = jce1, jce2
        zz1 = zeta(j,i)
        zlr = stdlrate(yeardayfrac(rcmtimer%idate),dayspy,mddom%xlat(j,i))
        ! zlr = -lrate
        tv = xtb(j,i,kz) * (d_one + ep1*xqb(j,i,kz)) + d_half * zz1 * zlr
        zz2 = egrav/(rgas*tv)
        p = xpsb(j,i) * exp(-zz1*zz2)
        xpaib(j,i) = (p/p00)**rovcp
      end do
    end do
  end subroutine paicompute

end module mod_atm_stub

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
