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
  use mod_memutil
  use mpi

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

  type diffx
    real(8) , pointer , dimension(:,:,:) :: difft
    real(8) , pointer , dimension(:,:,:) :: difuu
    real(8) , pointer , dimension(:,:,:) :: difuv
    real(8) , pointer , dimension(:,:,:) :: diffq
  end type diffx

  type v3dbound
    real(8) , pointer , dimension(:,:,:) :: b0
    real(8) , pointer , dimension(:,:,:) :: b1
    real(8) , pointer , dimension(:,:,:) :: nb , sb
    real(8) , pointer , dimension(:,:,:) :: eb , wb
    real(8) , pointer , dimension(:,:,:) :: nbt , sbt
    real(8) , pointer , dimension(:,:,:) :: ebt , wbt
  end type v3dbound

  type v2dbound
    real(8) , pointer , dimension(:,:) :: b0
    real(8) , pointer , dimension(:,:) :: b1
    real(8) , pointer , dimension(:,:) :: nb , sb
    real(8) , pointer , dimension(:,:) :: eb , wb
    real(8) , pointer , dimension(:,:) :: nbt , sbt
    real(8) , pointer , dimension(:,:) :: ebt , wbt
  end type v2dbound

  public :: atmstate , domain , surfpstate , surftstate , surfstate , slice
  public :: diffx , v2dbound , v3dbound

  type(domain) , public :: mddom
  type(atmstate) , public :: atm1 , atm2
  type(atmstate) , public :: atmx , atmc , aten , holtten , uwten
  type(surfpstate) , public :: sps1 , sps2
  type(surftstate) , public :: sts1 , sts2
  type(surfstate) , public :: sfsta
  type(slice) , public :: atms
  type(diffx) , public :: adf
  type(v3dbound) , public :: xtb , xqb , xub , xvb
  type(v2dbound) , public :: xpsb

  public :: allocate_mod_atm_interface , allocate_atmstate , allocate_domain
  public :: uvcross2dot

  real(8) , public , pointer , dimension(:,:) :: hgfact
  real(8) , public , pointer , dimension(:,:,:) :: sulfate
  real(8) , public , pointer, dimension(:,:,:) :: dstor
  real(8) , public , pointer, dimension(:,:,:) :: hstor
!
  real(8) , public , pointer , dimension(:,:) :: psc , pten , psd
  real(8) , public , pointer , dimension(:,:,:) :: phi , qdot , omega
!
  contains 
!
    subroutine allocate_v3dbound(xb,lband,ke,nsp)
      type(v3dbound) , intent(out) :: xb
      logical , intent(in) :: lband
      integer , intent(in) :: ke , nsp
      call getmem3d(xb%b0,1,iy,1,ke,1,jxp,'v3dbound:b0')
      call getmem3d(xb%b1,1,iy,1,ke,1,jxp,'v3dbound:b1')
      call getmem3d(xb%nb,1,nsp,1,ke,0,jxp+1,'v3dbound:nb')
      call getmem3d(xb%sb,1,nsp,1,ke,0,jxp+1,'v3dbound:sb')
      call getmem3d(xb%nbt,1,nsp,1,ke,0,jxp+1,'v3dbound:nbt')
      call getmem3d(xb%sbt,1,nsp,1,ke,0,jxp+1,'v3dbound:sbt')
      if ( .not. lband ) then
        call getmem3d(xb%eb,1,iy,1,ke,0,jxp+1,'v3dbound:eb')
        call getmem3d(xb%wb,1,iy,1,ke,0,jxp+1,'v3dbound:wb')
        call getmem3d(xb%ebt,1,iy,1,ke,0,jxp+1,'v3dbound:ebt')
        call getmem3d(xb%wbt,1,iy,1,ke,0,jxp+1,'v3dbound:wbt')
      end if
    end subroutine allocate_v3dbound
!
    subroutine allocate_v2dbound(xb,lband,nsp)
      type(v2dbound) , intent(out) :: xb
      logical , intent(in) :: lband
      integer , intent(in) :: nsp
      call getmem2d(xb%b0,1,iy,0,jxp+1,'v2dbound:b0')
      call getmem2d(xb%b1,1,iy,0,jxp+1,'v2dbound:b1')
      call getmem2d(xb%nb,1,nsp,0,jxp+1,'v2dbound:nb')
      call getmem2d(xb%sb,1,nsp,0,jxp+1,'v2dbound:sb')
      call getmem2d(xb%nbt,1,nsp,0,jxp+1,'v2dbound:nbt')
      call getmem2d(xb%sbt,1,nsp,0,jxp+1,'v2dbound:sbt')
      if ( .not. lband ) then
        call getmem2d(xb%eb,1,iy,0,jxp+1,'v2dbound:eb')
        call getmem2d(xb%wb,1,iy,0,jxp+1,'v2dbound:wb')
        call getmem2d(xb%ebt,1,iy,0,jxp+1,'v2dbound:ebt')
        call getmem2d(xb%wbt,1,iy,0,jxp+1,'v2dbound:wbt')
      end if
    end subroutine allocate_v2dbound
!
    subroutine allocate_atmstate(atm,ibltyp,lpar,ib,jb)
      implicit none
      logical , intent(in) :: lpar
      integer , intent(in) :: ibltyp
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
    subroutine allocate_diffx(dx)
      implicit none
      type(diffx) , intent(out) :: dx
      call getmem3d(dx%difft,1,iy,1,kz,1,jxp,'diffx:difft')
      call getmem3d(dx%difuu,1,iy,1,kz,1,jxp,'diffx:difuu')
      call getmem3d(dx%difuv,1,iy,1,kz,1,jxp,'diffx:difuv')
      call getmem3d(dx%diffq,1,iy,1,kz,1,jxp,'diffx:diffq')
    end subroutine allocate_diffx
!
    subroutine allocate_mod_atm_interface(lband,ibltyp)
!
      implicit none
      logical , intent(in) :: lband
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

      call allocate_surfpstate(sps1)
      call allocate_surfpstate(sps2)

      call allocate_surftstate(sts1)
      call allocate_surftstate(sts2)

      call allocate_surfstate(sfsta)

      call allocate_slice(atms)

      call allocate_diffx(adf)

      call allocate_v3dbound(xtb,lband,kz,nspgx)
      call allocate_v3dbound(xqb,lband,kz,nspgx)
      call allocate_v3dbound(xub,lband,kz,nspgd)
      call allocate_v3dbound(xvb,lband,kz,nspgd)
      call allocate_v2dbound(xpsb,lband,nspgx)

      if ( ehso4 ) then
        call getmem3d(sulfate,1,iy,1,kz,1,jxp,'mod_atm_interface:sulfate')
      end if

      call getmem3d(dstor,1,iy,0,jxp+1,1,nsplit,'mod_atm_interface:dstor')
      call getmem3d(hstor,1,iy,0,jxp+1,1,nsplit,'mod_atm_interface:hstor')
!
      call getmem2d(hgfact,1,iy,1,jxp,'mod_atm_interface:hgfact')
      call getmem3d(omega,1,iy,1,kz,1,jxp,'mod_atm_interface:omega')
      call getmem2d(psc,1,iy,1,jxp,'mod_atm_interface:psc')
      call getmem2d(pten,1,iy,1,jxp,'mod_atm_interface:pten')
      call getmem3d(phi,1,iy,1,kz,0,jxp,'mod_atm_interface:phi')
      call getmem2d(psd,1,iy,0,jxp+1,'mod_atm_interface:psd')
      call getmem3d(qdot,1,iy,1,kzp1,0,jxp+1,'mod_atm_interface:qdot')

    end subroutine allocate_mod_atm_interface 
!
! Takes an atmstate variable with u and v on the cross grid (the
! same grid as t, qv, qc, etc.) and interpolates the u and v to
! the dot grid.  This routine sheilds the user of the function
! from the need to worry about the details of the domain
! decomposition.  
!
! Written by Travis A. O'Brien 01/04/11.
!
! type(atmstate),intent(in) :: invar 
!                              An atmstate variable (see mod_atm_interface)
!                              that contains the u and v variables
!                              that need to be interpolated.  This
!                              variable is not modified by this
!                              routine.
!
! type(atmstate),intent(inout) :: outvar 
!                              An atmstate variable (see mod_atm_interface)
!                              that contains the u and v variables
!                              that will be overwritten by the
!                              interpolation of invar%u and invar%v.
!                              Only u and v are modified in this
!                              routine (t, qv, qc, and tke should
!                              remain unchanged).
!
    subroutine uvcross2dot(invar,outvar)
      implicit none
      type(atmstate) , intent(inout) :: invar
      type(atmstate) , intent(inout) :: outvar
      integer :: ib , ie , jb , je , i , j
      integer :: isendcount , ierr

      ! TODO:  It might make sense to encapsulate the following code
      ! in to a standard routine, since this boundary sending code is
      ! ubiquitous throughout the RegCM code and it is domain
      ! decomposition-dependent.

      ! Send the right-edge of the u/v tendencies to the left
      ! edge of the next process's u/v tendencies (so that
      ! invar%u(i,k,0) holds invar%u(i,k,jxp) of the parallel
      ! chunk next door)

      isendcount = iy*kz
      call mpi_sendrecv(invar%u(:,:,jxp),isendcount,mpi_real8,ieast,30, &
                        invar%u(:,:,0),isendcount,mpi_real8,iwest,30,   &
                        mycomm,mpi_status_ignore,ierr)
      call mpi_sendrecv(invar%v(:,:,jxp),isendcount,mpi_real8,ieast,31, &
                        invar%v(:,:,0),isendcount,mpi_real8,iwest,31,   &
                        mycomm,mpi_status_ignore,ierr)

      ! Set j-loop boundaries
      jb = jbegin
      je = jendx
      ! Set i-loop boundaries
      ib = 2
      ie = iym1

      !
      !     x     x     x     x     x     x
      !
      !        o     o     o     o     o 
      !         (i-1,j-1)     (i,j-1)            
      !     x     x     x-----x     x     x
      !                 |(i,j)|
      !        o     o  |  o  |  o     o 
      !                 |     |
      !     x     x     x-----x     x     x
      !           (i-1,j)     (i,j)
      !
      !        o     o     o     o     o 
      !
      !     x     x     x     x     x     x
      !

      ! Perform the bilinear interpolation necessary
      ! to put the u and v variables on the dot grid.

      do j = jb , je
        do i = ib , ie
          outvar%u(i,:,j) =  outvar%u(i,:,j) +             &
            d_rfour*(invar%u(i,:,j) + invar%u(i,:,j-1) +   &
                     invar%u(i-1,:,j) + invar%u(i-1,:,j-1))
          outvar%v(i,:,j) =  outvar%v(i,:,j) +             &
            d_rfour*(invar%v(i,:,j) + invar%v(i,:,j-1) +   &
                     invar%v(i-1,:,j) + invar%v(i-1,:,j-1))
        end do
      end do
    end subroutine uvcross2dot
!
end module mod_atm_interface
