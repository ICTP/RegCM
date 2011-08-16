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
  use mod_runparams
  use mod_memutil

  private

!
! Storage for all the 3d prognostic variables in two
!     timesteps and all the 2d variables and constants
!
  type domain
    real(8) , pointer , dimension(:,:) :: ht
    real(8) , pointer , dimension(:,:) :: lndcat
    real(8) , pointer , dimension(:,:) :: xlat
    real(8) , pointer , dimension(:,:) :: xlon
    real(8) , pointer , dimension(:,:) :: msfx
    real(8) , pointer , dimension(:,:) :: msfd
    real(8) , pointer , dimension(:,:) :: coriol
  end type domain

  type atmstate
    real(8) , pointer , dimension(:,:,:) :: u
    real(8) , pointer , dimension(:,:,:) :: v
    real(8) , pointer , dimension(:,:,:) :: t
    real(8) , pointer , dimension(:,:,:) :: qv
    real(8) , pointer , dimension(:,:,:) :: qc
    real(8) , pointer , dimension(:,:,:) :: tke
  end type atmstate

  type surfpstate
    real(8) , pointer , dimension(:,:) :: ps
    real(8) , pointer , dimension(:,:) :: pdot
  end type surfpstate

  type surftstate
    real(8) , pointer , dimension(:,:) :: tg
  end type surftstate

  type surfstate
    real(8) , pointer , dimension(:,:) :: rainc
    real(8) , pointer , dimension(:,:) :: rainnc
    real(8) , pointer , dimension(:,:) :: hfx
    real(8) , pointer , dimension(:,:) :: qfx
    real(8) , pointer , dimension(:,:) :: tgbb
    real(8) , pointer , dimension(:,:) :: zpbl
    real(8) , pointer , dimension(:,:) :: uvdrag
  end type surfstate

  type slice
    real(8) , pointer , dimension(:,:,:) :: tb3d
    real(8) , pointer , dimension(:,:,:) :: thx3d
    real(8) , pointer , dimension(:,:,:) :: pb3d
    real(8) , pointer , dimension(:,:,:) :: rhob3d
    real(8) , pointer , dimension(:,:,:) :: ubx3d
    real(8) , pointer , dimension(:,:,:) :: vbx3d
    real(8) , pointer , dimension(:,:,:) :: ubd3d
    real(8) , pointer , dimension(:,:,:) :: vbd3d
    real(8) , pointer , dimension(:,:,:) :: rhb3d
    real(8) , pointer , dimension(:,:,:) :: qvb3d
    real(8) , pointer , dimension(:,:,:) :: qsb3d
    real(8) , pointer , dimension(:,:,:) :: qcb3d
    real(8) , pointer , dimension(:,:,:,:) :: chib3d
  end type slice

  type(domain) , public :: mddom
  type(atmstate) , public :: atm1 , atm2
  type(atmstate) , public :: atmx , atmc , aten , holtten
  type(surfpstate) , public :: sps1 , sps2
  type(surftstate) , public :: sts1 , sts2
  type(surfstate) , public :: sfsta
  type(slice) , public :: atms

  public :: atmstate , domain , surfpstate , surftstate , surfstate , slice
  public :: allocate_mod_atm_interface , allocate_atmstate , allocate_domain

  real(8) , public , pointer , dimension(:,:) :: hgfact
  integer , public , pointer , dimension(:,:) :: cucontrol
  real(8) , public , pointer , dimension(:,:,:) :: sulfate
  real(8) , public , pointer, dimension(:,:,:) :: dstor
  real(8) , public , pointer, dimension(:,:,:) :: hstor
!
  real(8) , public , pointer , dimension(:,:,:) :: diffq , difft , &
                                                   difuu , difuv
  real(8) , public , pointer , dimension(:,:) :: psc , pten , psd
  real(8) , public , pointer , dimension(:,:,:) :: phi , qdot , omega
!
  contains 
!
    subroutine allocate_atmstate(atm,lpar,ib,jb)
      implicit none
      logical , intent(in) :: lpar
      integer , intent(in) :: ib , jb
      type(atmstate) , intent(out) :: atm
      integer :: is , ie , js , je
      if (lpar) then
        is = 1-ib
        ie = iy+ib
        js = 1-jb
        je = jxp+jb
        call getmem3d(atm%u,is,ie,1,kz,js,je,'atmstate:u')
        call getmem3d(atm%v,is,ie,1,kz,js,je,'atmstate:v')
        call getmem3d(atm%t,is,ie,1,kz,js,je,'atmstate:t')
        call getmem3d(atm%qv,is,ie,1,kz,js,je,'atmstate:qv')
        call getmem3d(atm%qc,is,ie,1,kz,js,je,'atmstate:qc')
        if ( ibltyp == 2 .or. ibltyp == 99 ) then
          call getmem3d(atm%tke,is,ie,1,kzp1,js,je,'atmstate:tke')
        end if
      else
        call getmem3d(atm%u,1,iy,1,kz,1,jx,'atmstate:u')
        call getmem3d(atm%v,1,iy,1,kz,1,jx,'atmstate:v')
        call getmem3d(atm%t,1,iy,1,kz,1,jx,'atmstate:t')
        call getmem3d(atm%qv,1,iy,1,kz,1,jx,'atmstate:qv')
        call getmem3d(atm%qc,1,iy,1,kz,1,jx,'atmstate:qc')
        if ( ibltyp == 2 .or. ibltyp == 99 ) then
          call getmem3d(atm%tke,1,iy,1,kzp1,1,jx,'atmstate:tke')
        end if
      end if
    end subroutine allocate_atmstate
!
    subroutine allocate_surfpstate(sps)
      implicit none
      type(surfpstate) , intent(out) :: sps
      call getmem2d(sps%ps,1,iy,-1,jxp+2,'surfpstate:ps')
      call getmem2d(sps%pdot,1,iy,-1,jxp+2,'surfpstate:pdot')
    end subroutine allocate_surfpstate
!
    subroutine allocate_surftstate(sts)
      implicit none
      type(surftstate) , intent(out) :: sts
      call getmem2d(sts%tg,1,iy,1,jxp+1,'surftstate:tg')
    end subroutine allocate_surftstate
!
    subroutine allocate_domain(dom,lpar)
      implicit none
      logical , intent(in) :: lpar
      type(domain) , intent(out) :: dom

      if (lpar) then
        call getmem2d(dom%ht,1,iy,0,jxp+1,'mod_atm_interface:ht')
        call getmem2d(dom%lndcat,1,iy,1,jxp,'mod_atm_interface:lndcat')
        call getmem2d(dom%xlat,1,iy,1,jxp,'mod_atm_interface:xlat')
        call getmem2d(dom%xlon,1,iy,1,jxp,'mod_atm_interface:xlon')
        call getmem2d(dom%msfx,1,iy,-1,jxp+2,'mod_atm_interface:msfx')
        call getmem2d(dom%msfd,1,iy,-1,jxp+2,'mod_atm_interface:msfd')
        call getmem2d(dom%coriol,1,iy,1,jxp,'mod_atm_interface:f')
      else
        call getmem2d(dom%ht,1,iy,1,jx,'mod_atm_interface:ht')
        call getmem2d(dom%lndcat,1,iy,1,jx,'mod_atm_interface:lndcat')
        call getmem2d(dom%xlat,1,iy,1,jx,'mod_atm_interface:xlat')
        call getmem2d(dom%xlon,1,iy,1,jx,'mod_atm_interface:xlon')
        call getmem2d(dom%msfx,1,iy,1,jx,'mod_atm_interface:msfx')
        call getmem2d(dom%msfd,1,iy,1,jx,'mod_atm_interface:msfd')
        call getmem2d(dom%coriol,1,iy,1,jx,'mod_atm_interface:f')
      end if
    end subroutine allocate_domain
!
    subroutine allocate_surfstate(sfs)
      implicit none
      type(surfstate) , intent(out) :: sfs
      call getmem2d(sfs%hfx,1,iy,1,jxp,'surfstate:hfx')
      call getmem2d(sfs%qfx,1,iy,1,jxp,'surfstate:qfx')
      call getmem2d(sfs%rainc,1,iy,1,jxp,'surfstate:rainc')
      call getmem2d(sfs%rainnc,1,iy,1,jxp,'surfstate:rainnc')
      call getmem2d(sfs%tgbb,1,iy,1,jxp,'surfstate:tgbb')
      call getmem2d(sfs%zpbl,1,iy,1,jxp,'surfstate:zpbl')
      call getmem2d(sfs%uvdrag,1,iy,0,jxp,'surfstate:uvdrag')
    end subroutine allocate_surfstate
!
    subroutine allocate_slice(ax)
      implicit none
      type(slice) , intent(out) :: ax
      call getmem3d(ax%pb3d,1,iy,1,kz,1,jxp,'slice:pb3d')
      call getmem3d(ax%qsb3d,1,iy,1,kz,1,jxp,'slice:qsb3d')
      call getmem3d(ax%rhb3d,1,iy,1,kz,1,jxp,'slice:rhb3d')
      call getmem3d(ax%rhob3d,1,iy,1,kz,1,jxp,'slice:rhob3d')
      call getmem3d(ax%ubx3d,1,iy,1,kz,1,jxp,'slice:ubx3d')
      call getmem3d(ax%vbx3d,1,iy,1,kz,1,jxp,'slice:vbx3d')
      call getmem3d(ax%thx3d,1,iy,1,kz,1,jxp,'slice:thx3d')
      call getmem3d(ax%qcb3d,1,iy,1,kz,-1,jxp+2,'slice:qcb3d')
      call getmem3d(ax%qvb3d,1,iy,1,kz,-1,jxp+2,'slice:qvb3d')
      call getmem3d(ax%tb3d,1,iy,1,kz,-1,jxp+2,'slice:tb3d')
      call getmem3d(ax%ubd3d,1,iy,1,kz,-1,jxp+2,'slice:ubd3d')
      call getmem3d(ax%vbd3d,1,iy,1,kz,-1,jxp+2,'slice:vbd3d')
      if ( ichem == 1 ) then
        call getmem4d(ax%chib3d,1,iy,1,kz,-1,jxp+2,1,ntr,'slice:chib3d')
      end if
    end subroutine allocate_slice
!
    subroutine allocate_mod_atm_interface
      implicit none

      call allocate_domain(mddom,.true.)

      call allocate_atmstate(atm1,.true.,0,2)
      call allocate_atmstate(atm2,.true.,0,2)
      call allocate_atmstate(atmx,.true.,0,1)
      call allocate_atmstate(atmc,.true.,0,0)
      call allocate_atmstate(aten,.true.,0,0)
      if ( ibltyp == 99 ) then
        call allocate_atmstate(holtten,.true.,0,0)
      end if

      call allocate_surfpstate(sps1)
      call allocate_surfpstate(sps2)

      call allocate_surftstate(sts1)
      call allocate_surftstate(sts2)

      call allocate_surfstate(sfsta)

      call allocate_slice(atms)

      if (icup == 99 .or. icup == 98) then
        call getmem2d(cucontrol,1,iy,1,jxp,'mod_atm_interface:cucontrol')
      end if

      if ( ehso4 ) then
        call getmem3d(sulfate,1,iy,1,kz,1,jxp,'mod_atm_interface:sulfate')
      end if

      call getmem3d(dstor,1,iy,0,jxp+1,1,nsplit,'mod_atm_interface:dstor')
      call getmem3d(hstor,1,iy,0,jxp+1,1,nsplit,'mod_atm_interface:hstor')
!
      call getmem2d(hgfact,1,iy,1,jxp,'mod_atm_interface:hgfact')
      call getmem3d(diffq,1,iy,1,kz,1,jxp,'mod_atm_interface:diffq')
      call getmem3d(difft,1,iy,1,kz,1,jxp,'mod_atm_interface:difft')
      call getmem3d(difuu,1,iy,1,kz,1,jxp,'mod_atm_interface:difuu')
      call getmem3d(difuv,1,iy,1,kz,1,jxp,'mod_atm_interface:difuv')
      call getmem3d(omega,1,iy,1,kz,1,jxp,'mod_atm_interface:omega')
      call getmem2d(psc,1,iy,1,jxp,'mod_atm_interface:psc')
      call getmem2d(pten,1,iy,1,jxp,'mod_atm_interface:pten')
      call getmem3d(phi,1,iy,1,kz,0,jxp,'mod_atm_interface:phi')
      call getmem2d(psd,1,iy,0,jxp+1,'mod_atm_interface:psd')
      call getmem3d(qdot,1,iy,1,kzp1,0,jxp+1,'mod_atm_interface:qdot')

    end subroutine allocate_mod_atm_interface 
!
end module mod_atm_interface
