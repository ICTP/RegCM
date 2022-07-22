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

  implicit none

  private

  type(domain) , public :: mddom
  type(domain_subgrid) , public :: mdsub
  type(surfstate) , public :: sfs

  type(v3dbound) , public :: xtb , xqb , xub , xvb , xlb , xib , xppb , xwwb
  type(v2dbound) , public :: xpsb , xtsb , xpaib

  type(lm_exchange) , public :: lm
  type(lm_state) , public :: lms

  integer(ik4) :: ix1 , ix2 , jx1 , jx2
  integer(ik4) :: id1 , id2 , jd1 , jd2

  real(rkx) , public :: rdnnsg
  real(rkx) :: rdtbdy

  logical , parameter :: cross = .false.
  logical , parameter :: dot = .true.

  real(rkx) , dimension(:,:) , pointer , public :: patm
  real(rkx) , dimension(:,:) , pointer , public :: tatm
  real(rkx) , dimension(:,:) , pointer , public :: uatm
  real(rkx) , dimension(:,:) , pointer , public :: vatm
  real(rkx) , dimension(:,:) , pointer , public :: thatm
  real(rkx) , dimension(:,:) , pointer , public :: qvatm
  real(rkx) , dimension(:,:) , pointer , public :: zatm
  real(rkx) , dimension(:,:) , pointer , public :: rho
  real(rkx) , dimension(:,:) , pointer , public :: ps
  real(rkx) , dimension(:,:) , pointer , public :: tp
  real(rkx) , dimension(:,:) , pointer , public :: coszrs
  real(rkx) , dimension(:,:) , pointer , public :: fsw
  real(rkx) , dimension(:,:) , pointer , public :: flw
  real(rkx) , dimension(:,:) , pointer , public :: flwd
  real(rkx) , dimension(:,:) , pointer , public :: solar
  real(rkx) , dimension(:,:) , pointer , public :: pptc
  real(rkx) , dimension(:,:) , pointer , public :: pptnc
  real(rkx) , dimension(:,:) , pointer , public :: swdir
  real(rkx) , dimension(:,:) , pointer , public :: swdif
  real(rkx) , dimension(:,:) , pointer , public :: lwdir
  real(rkx) , dimension(:,:) , pointer , public :: lwdif
  real(rkx) , dimension(:,:) , pointer , public :: totc
  real(rkx) , dimension(:,:) , pointer , public :: zeta , fmzf

  public :: allocate_mod_atm_interface , allocate_surface_model
  public :: setup_model_indexes
  public :: init_bdy , bdyin , atmval

  interface timeint
    module procedure timeint2 , timeint3
  end interface timeint

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
        call getmem3d(xb%b0,jce1,jce2,ice1,ice2,1,ke,'v3dbound:b0')
        call getmem3d(xb%b1,jce1,jce2,ice1,ice2,1,ke,'v3dbound:b1')
        call getmem3d(xb%bt,jce1,jce2,ice1,ice2,1,ke,'v3dbound:bt')
      end if
    end subroutine allocate_v3dbound

    subroutine allocate_v2dbound(xb,ldot)
      implicit none
      type(v2dbound) , intent(out) :: xb
      logical , intent(in) :: ldot
      if ( ldot ) then
        call getmem2d(xb%b0,jde1,jde2,ide1,ide2,'v2dbound:b0')
        call getmem2d(xb%b1,jde1,jde2,ide1,ide2,'v2dbound:b1')
        call getmem2d(xb%bt,jde1,jde2,ide1,ide2,'v2dbound:bt')
      else
        call getmem2d(xb%b0,jce1,jce2,ice1,ice2,'v2dbound:b0')
        call getmem2d(xb%b1,jce1,jce2,ice1,ice2,'v2dbound:b1')
        call getmem2d(xb%bt,jce1,jce2,ice1,ice2,'v2dbound:bt')
      end if
    end subroutine allocate_v2dbound

    subroutine allocate_domain(dom)
      implicit none
      type(domain) , intent(out) :: dom
      call getmem2d(dom%ht,jde1,jde2,ide1,ide2,'storage:ht')
      call getmem2d(dom%lndcat,jde1,jde2,ide1,ide2,'storage:lndcat')
      call getmem2d(dom%lndtex,jde1,jde2,ide1,ide2,'storage:lndtex')
      call getmem2d(dom%xlat,jde1,jde2,ide1,ide2,'storage:xlat')
      call getmem2d(dom%xlon,jde1,jde2,ide1,ide2,'storage:xlon')
      call getmem2d(dom%mask,jde1,jde2,ide1,ide2,'storage:mask')
      if ( idynamic == 3 ) then
        call getmem2d(dom%msfx,jde1,jde2,ide1,ide2,'storage:msfx')
        call getmem2d(dom%msfu,jde1,jde2,ide1,ide2,'storage:msfu')
        call getmem2d(dom%msfv,jde1,jde2,ide1,ide2,'storage:msfv')
        call getmem2d(dom%hx,jde1,jde2,ice1,ice2,'storage:hx')
        call getmem2d(dom%hy,jce1,jce2,ide1,ide2,'storage:hy')
        call getmem2d(dom%ulat,jde1,jde2,ide1,ide2,'storage:ulat')
        call getmem2d(dom%ulon,jde1,jde2,ide1,ide2,'storage:ulon')
        call getmem2d(dom%vlat,jde1,jde2,ide1,ide2,'storage:vlat')
        call getmem2d(dom%vlon,jde1,jde2,ide1,ide2,'storage:vlon')
        call getmem2d(dom%coriol,jde1,jde2,ide1,ide2,'storage:f')
      else
        call getmem2d(dom%dlat,jde1,jde2,ide1,ide2,'storage:dlat')
        call getmem2d(dom%dlon,jde1,jde2,ide1,ide2,'storage:dlon')
        call getmem2d(dom%msfx,jd1,jd2,id1,id2,'storage:msfx')
        call getmem2d(dom%msfd,jd1,jd2,id1,id2,'storage:msfd')
        call getmem2d(dom%coriol,jde1,jde2,ide1,ide2,'storage:f')
      end if
      call getmem2d(dom%snowam,jde1,jde2,ide1,ide2,'storage:snowam')
      call getmem2d(dom%smoist,jde1,jde2,ide1,ide2,'storage:smoist')
      call getmem3d(dom%rmoist,jde1,jde2,ide1,ide2, &
                    1,num_soil_layers,'storage:rmoist')
      call getmem2d(dom%ldmsk,jci1,jci2,ici1,ici2,'storage:ldmsk')
      call getmem2d(dom%iveg,jci1,jci2,ici1,ici2,'storage:iveg')
      call getmem2d(dom%itex,jci1,jci2,ici1,ici2,'storage:itex')
      call getmem2d(dom%xmsf,jdi1,jdi2,idi1,idi2,'storage:xmsf')
      call getmem2d(dom%dmsf,jdi1,jdi2,idi1,idi2,'storage:dmsf')
    end subroutine allocate_domain

    subroutine allocate_domain_subgrid(sub)
      implicit none
      type(domain_subgrid) , intent(out) :: sub
      call getmem3d(sub%ht,1,nnsg,jde1,jde2,ide1,ide2,'storage:ht')
      call getmem3d(sub%lndcat,1,nnsg,jde1,jde2,ide1,ide2,'storage:lndcat')
      call getmem3d(sub%lndtex,1,nnsg,jde1,jde2,ide1,ide2,'storage:lndtex')
      call getmem3d(sub%xlat,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlat')
      call getmem3d(sub%xlon,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlon')
      call getmem3d(sub%mask,1,nnsg,jde1,jde2,ide1,ide2,'storage:xlon')
      call getmem3d(sub%ldmsk,1,nnsg,jci1,jci2,ici1,ici2,'storage:ldmsk')
      call getmem3d(sub%iveg,1,nnsg,jci1,jci2,ici1,ici2,'storage:iveg')
      call getmem3d(sub%itex,1,nnsg,jci1,jci2,ici1,ici2,'storage:itex')
    end subroutine allocate_domain_subgrid

    subroutine allocate_surfstate(sfs)
      implicit none
      type(surfstate) , intent(out) :: sfs
      call getmem2d(sfs%psa,jce1,jce2,ice1,ice2,'surf:psa')
      call getmem2d(sfs%psdota,jde1,jde2,ide1,ide2,'surf:psdota')
      call getmem2d(sfs%psb,jx1,jx2,ix1,ix2,'surf:psb')
      call getmem2d(sfs%psdotb,jd1,jd2,id1,id2,'surf:psdotb')
      call getmem2d(sfs%psc,jce1,jce2,ice1,ice2,'surf:psc')
      call getmem2d(sfs%tg,jci1,jci2,ici1,ici2,'surf:tg')
      call getmem2d(sfs%hfx,jci1,jci2,ici1,ici2,'surf:hfx')
      call getmem2d(sfs%qfx,jci1,jci2,ici1,ici2,'surf:qfx')
      call getmem2d(sfs%rainc,jci1,jci2,ici1,ici2,'surf:rainc')
      call getmem2d(sfs%rainnc,jci1,jci2,ici1,ici2,'surf:rainnc')
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
    end subroutine allocate_surfstate

    subroutine allocate_mod_atm_interface
      implicit none

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
        call getmem2d(zeta,jce1,jce2,ice1,ice2,'lm:zeta')
        call getmem2d(fmzf,jce1,jce2,ice1,ice2,'lm:fmzf')
      end if
    end subroutine allocate_mod_atm_interface

    subroutine allocate_surface_model
      implicit none

      rdnnsg = d_one/real(nnsg,rkx)

      call getmem3d(lms%sent,1,nnsg,jci1,jci2,ici1,ici2,'lm:sent')
      call getmem3d(lms%evpr,1,nnsg,jci1,jci2,ici1,ici2,'lm:evpr')
      call getmem3d(lms%deltat,1,nnsg,jci1,jci2,ici1,ici2,'lm:deltat')
      call getmem3d(lms%deltaq,1,nnsg,jci1,jci2,ici1,ici2,'lm:deltaq')
      call getmem3d(lms%drag,1,nnsg,jci1,jci2,ici1,ici2,'lm:drag')
      call getmem3d(lms%ustar,1,nnsg,jci1,jci2,ici1,ici2,'lm:ustar')
      call getmem3d(lms%w10m,1,nnsg,jci1,jci2,ici1,ici2,'lm:w10m')
      call getmem3d(lms%zo,1,nnsg,jci1,jci2,ici1,ici2,'lm:zo')
      call getmem3d(lms%rhoa,1,nnsg,jci1,jci2,ici1,ici2,'lm:rho')
      call getmem3d(lms%lncl,1,nnsg,jci1,jci2,ici1,ici2,'lm:lncl')
      call getmem3d(lms%prcp,1,nnsg,jci1,jci2,ici1,ici2,'lm:prcp')
      call getmem3d(lms%snwm,1,nnsg,jci1,jci2,ici1,ici2,'lm:snwm')
      call getmem3d(lms%trnof,1,nnsg,jci1,jci2,ici1,ici2,'lm:trnof')
      call getmem3d(lms%srnof,1,nnsg,jci1,jci2,ici1,ici2,'lm:srnof')
      call getmem3d(lms%xlai,1,nnsg,jci1,jci2,ici1,ici2,'lm:xlai')
      call getmem3d(lms%sfcp,1,nnsg,jci1,jci2,ici1,ici2,'lm:sfcp')
      call getmem3d(lms%q2m,1,nnsg,jci1,jci2,ici1,ici2,'lm:q2m')
      call getmem3d(lms%t2m,1,nnsg,jci1,jci2,ici1,ici2,'lm:t2m')
      call getmem3d(lms%u10m,1,nnsg,jci1,jci2,ici1,ici2,'lm:u10m')
      call getmem3d(lms%v10m,1,nnsg,jci1,jci2,ici1,ici2,'lm:v10m')
      call getmem3d(lms%ram1,1,nnsg,jci1,jci2,ici1,ici2,'lm:ram1')
      call getmem3d(lms%rah1,1,nnsg,jci1,jci2,ici1,ici2,'lm:rah1')
      call getmem3d(lms%br,1,nnsg,jci1,jci2,ici1,ici2,'lm:br')
      call getmem3d(lms%taux,1,nnsg,jci1,jci2,ici1,ici2,'lm:taux')
      call getmem3d(lms%tauy,1,nnsg,jci1,jci2,ici1,ici2,'lm:tauy')
      call getmem3d(lms%wt,1,nnsg,jci1,jci2,ici1,ici2,'lm:wt')
      call getmem3d(lms%swalb,1,nnsg,jci1,jci2,ici1,ici2,'lm:swalb')
      call getmem3d(lms%lwalb,1,nnsg,jci1,jci2,ici1,ici2,'lm:lwalb')
      call getmem3d(lms%swdiralb,1,nnsg,jci1,jci2,ici1,ici2,'lm:swdiralb')
      call getmem3d(lms%lwdiralb,1,nnsg,jci1,jci2,ici1,ici2,'lm:lwdiralb')
      call getmem3d(lms%swdifalb,1,nnsg,jci1,jci2,ici1,ici2,'lm:swdifalb')
      call getmem3d(lms%lwdifalb,1,nnsg,jci1,jci2,ici1,ici2,'lm:lwdifalb')

      call getmem3d(lms%gwet,1,nnsg,jci1,jci2,ici1,ici2,'lm:gwet')
      call getmem4d(lms%sw,1,nnsg,jci1,jci2,ici1,ici2,1,num_soil_layers,'lm:sw')

      call assignpnt(lms%sw,lms%ssw,1)
      call assignpnt(lms%sw,lms%rsw,2)
      call assignpnt(lms%sw,lms%tsw,3)
      call getmem3d(lms%tgbb,1,nnsg,jci1,jci2,ici1,ici2,'lm:tgbb')
      call getmem3d(lms%tgrd,1,nnsg,jci1,jci2,ici1,ici2,'lm:tgrd')
      call getmem3d(lms%tgbrd,1,nnsg,jci1,jci2,ici1,ici2,'lm:tgbrd')
      call getmem3d(lms%tlef,1,nnsg,jci1,jci2,ici1,ici2,'lm:tlef')
      call getmem3d(lms%taf,1,nnsg,jci1,jci2,ici1,ici2,'lm:taf')
      call getmem3d(lms%sigf,1,nnsg,jci1,jci2,ici1,ici2,'lm:sigf')
      call getmem3d(lms%sfice,1,nnsg,jci1,jci2,ici1,ici2,'lm:sfice')
      call getmem3d(lms%snag,1,nnsg,jci1,jci2,ici1,ici2,'lm:snag')
      call getmem3d(lms%ldew,1,nnsg,jci1,jci2,ici1,ici2,'lm:ldew')
      call getmem3d(lms%sncv,1,nnsg,jci1,jci2,ici1,ici2,'lm:sncv')
      call getmem3d(lms%scvk,1,nnsg,jci1,jci2,ici1,ici2,'lm:scvk')
      call getmem3d(lms%um10,1,nnsg,jci1,jci2,ici1,ici2,'lm:um10')
      call getmem3d(lms%emisv,1,nnsg,jci1,jci2,ici1,ici2,'lm:emisv')
      call getmem4d(lms%vocemiss,1,nnsg,jci1,jci2,ici1,ici2,1,ntr,'lm:vocemiss')
      call getmem4d(lms%dustemiss,1,nnsg,jci1,jci2,ici1,ici2,1,4,'lm:dustemiss')
      call getmem4d(lms%sw_vol,1,nnsg,jci1,jci2, &
                                      ici1,ici2,1,num_soil_layers,'lm:sw_vol')
      call getmem4d(lms%tsoi,1,nnsg,jci1,jci2, &
                                    ici1,ici2,1,num_soil_layers,'lm:tsoi')

      call getmem2d(tatm,jci1,jci2,ici1,ici2,'lm:tatm')
      call getmem2d(patm,jci1,jci2,ici1,ici2,'lm:patm')
      call getmem2d(uatm,jci1,jci2,ici1,ici2,'lm:uatm')
      call getmem2d(vatm,jci1,jci2,ici1,ici2,'lm:vatm')
      call getmem2d(thatm,jci1,jci2,ici1,ici2,'lm:thatm')
      call getmem2d(qvatm,jci1,jci2,ici1,ici2,'lm:qvatm')
      call getmem2d(zatm,jci1,jci2,ici1,ici2,'lm:zatm')
      call getmem2d(rho,jci1,jci2,ici1,ici2,'lm:rho')
      call getmem2d(ps,jci1,jci2,ici1,ici2,'lm:ps')
      call getmem2d(tp,jci1,jci2,ici1,ici2,'lm:tp')
      call getmem2d(coszrs,jci1,jci2,ici1,ici2,'lm:coszrs')
      call getmem2d(flw,jci1,jci2,ici1,ici2,'lm:flw')
      call getmem2d(fsw,jci1,jci2,ici1,ici2,'lm:fsw')
      call getmem2d(flwd,jci1,jci2,ici1,ici2,'lm:flwd')
      call getmem2d(solar,jci1,jci2,ici1,ici2,'lm:solar')
      call getmem2d(pptc,jci1,jci2,ici1,ici2,'lm:pptc')
      call getmem2d(pptnc,jci1,jci2,ici1,ici2,'lm:pptnc')
      call getmem2d(swdir,jci1,jci2,ici1,ici2,'lm:swdir')
      call getmem2d(swdif,jci1,jci2,ici1,ici2,'lm:swdif')
      call getmem2d(lwdir,jci1,jci2,ici1,ici2,'lm:lwdir')
      call getmem2d(lwdif,jci1,jci2,ici1,ici2,'lm:lwdif')
      call getmem2d(totc,jci1,jci2,ici1,ici2,'lm:totc')
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
      implicit none
      character(len=32) :: appdat
      type (rcm_time_and_date) :: icbc_date
      integer(ik4) :: i , j , datefound
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'init_bdy'
      integer(ik4) , save :: idindx = 0
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
        write (stdout,*) 'READY  BC from     ' , &
              tochar10(bdydate1) , ' to ' , tochar10(bdydate2)
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
        do i = ici1 , ici2
          do j = jci1 , jci2
            sfs%tg(j,i) = xtsb%b0(j,i)
          end do
        end do
      end if

#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine init_bdy

    subroutine bdyin
      implicit none
      integer(ik4) :: datefound
      character(len=32) :: appdat
#ifdef DEBUG
      character(len=dbgslen) :: subroutine_name = 'bdyin'
      integer(ik4) , save :: idindx = 0
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
        write (stdout,*) 'READY  BC from     ' , &
              tochar10(bdydate1) , ' to ' , tochar10(bdydate2)
      end if

      bdydate1 = bdydate2

#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
    end subroutine bdyin

    subroutine atmval
      implicit none
      integer(ik4) :: i , j
      real(rkx) :: psb , cell , zq , xt

      xt = xbctime

      do i = ici1 , ici2
        do j = jci1 , jci2
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
      implicit none
      real(rkx) , pointer , dimension(:,:) , intent(in) :: a , b
      real(rkx) , pointer , dimension(:,:) , intent(inout) :: c
      integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
      integer(ik4) :: i , j
      do i = i1 , i2
        do j = j1 , j2
          c(j,i) = (a(j,i)-b(j,i))*rdtbdy
        end do
      end do
    end subroutine timeint2

    subroutine timeint3(a,b,c,j1,j2,i1,i2,k1,k2)
      implicit none
      real(rkx) , pointer , dimension(:,:,:) , intent(in) :: a , b
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: c
      integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
      integer(ik4) :: i , j , k
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            c(j,i,k) = (a(j,i,k)-b(j,i,k))*rdtbdy
          end do
        end do
      end do
    end subroutine timeint3

  subroutine paicompute(xpsb,xtb,xqb,xpaib)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xpsb
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: xtb , xqb
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: xpaib
    real(rkx) :: p , tv , zz1 , zz2 , zlr
    integer(ik4) :: i , j
    ! Hydrostatic initialization of pai
    do i = ice1 , ice2
      do j = jce1 , jce2
        zz1 = zeta(j,i)
        zlr = stdlrate(yeardayfrac(rcmtimer%idate),mddom%xlat(j,i))
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
